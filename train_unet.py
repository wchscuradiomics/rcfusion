import os
import random
import cv2
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import segmentation_models_pytorch as smp

def set_seed(seed):
  random.seed(seed)
  np.random.seed(seed)

  torch.manual_seed(seed)
  torch.cuda.manual_seed(seed)
  torch.cuda.manual_seed_all(seed)

  torch.backends.cudnn.deterministic = True
  torch.backends.cudnn.benchmark = False

# =========================
# 1. Dataset（带正负标记）
# =========================
class BrainCTDataset(Dataset):
  def __init__(self, img_dir, mask_dir):
    self.img_files = sorted([f for f in os.listdir(img_dir) if f.endswith(".png")])
    self.mask_files = sorted([f for f in os.listdir(mask_dir) if f.endswith(".png")])
    self.img_dir = img_dir
    self.mask_dir = mask_dir

    # 预计算是否含出血
    self.has_lesion = []
    for m in self.mask_files:
      mask = cv2.imread(os.path.join(mask_dir, m), 0)
      self.has_lesion.append(1 if np.any(mask == 1) else 0)

  def __len__(self):
    return len(self.img_files)

  def __getitem__(self, idx):
    img = cv2.imread(os.path.join(self.img_dir, self.img_files[idx]), 0).astype(np.float32) / 255.0
    mask = cv2.imread(os.path.join(self.mask_dir, self.mask_files[idx]), 0).astype(np.float32)
    mask = (mask == 1).astype(np.float32)

    img = torch.tensor(img[np.newaxis, :, :], dtype=torch.float32)
    mask = torch.tensor(mask[np.newaxis, :, :], dtype=torch.float32)

    return img, mask, self.has_lesion[idx]

# =========================
# 2. Loss（稳定版）
# =========================
class FocalTverskyLoss(nn.Module):
  def __init__(self, alpha=0.85, beta=0.15, gamma=0.75, smooth=1e-6):
    super().__init__()
    self.alpha = alpha
    self.beta = beta
    self.gamma = gamma
    self.smooth = smooth

  def forward(self, pred, target):
    pred = torch.sigmoid(pred)

    tp = (pred * target).sum(dim=(2,3))
    fp = ((1 - target) * pred).sum(dim=(2,3))
    fn = (target * (1 - pred)).sum(dim=(2,3))

    tversky = (tp + self.smooth) / (tp + self.alpha*fp + self.beta*fn + self.smooth)
    return ((1 - tversky) ** self.gamma).mean()

# =========================
# 3. Metrics
# =========================
def dice_iou(pred, mask):
  inter = (pred * mask).sum()
  dice = 2 * inter / (pred.sum() + mask.sum() + 1e-6)
  iou = inter / (pred.sum() + mask.sum() - inter + 1e-6)
  return dice.item(), iou.item()

# =========================
# 4. Evaluation（核心修复）
# =========================
def evaluate(model, loader, device, threshold=0.5, min_area=25):
  model.eval()

  all_d, all_i = [], []
  pos_d, pos_i = [], []

  fp_count = 0
  neg_count = 0

  with torch.no_grad():

    for img, mask, has_lesion in loader:

      img = img.to(device)
      mask = mask.to(device)

      pred = torch.sigmoid(model(img))

      pred_np = pred.squeeze().cpu().numpy()

      pred_bin = (pred_np > threshold).astype(np.uint8)

      num_labels, labels, stats, _ = \
        cv2.connectedComponentsWithStats(
          pred_bin,
          connectivity=8
        )

      for i in range(1, num_labels):

        area = stats[i, cv2.CC_STAT_AREA]

        if area < min_area:
          pred_bin[labels == i] = 0

      pred_bin = torch.tensor(
        pred_bin,
        dtype=torch.float32,
        device=device
      ).unsqueeze(0).unsqueeze(0)

      d, iou = dice_iou(pred_bin, mask)

      all_d.append(d)
      all_i.append(iou)

      if mask.sum() > 0:
        pos_d.append(d)
        pos_i.append(iou)

      if mask.sum() == 0:
        neg_count += 1

        if pred_bin.sum() > 50:
          fp_count += 1

  print("==============================")
  print("Evaluation Result")
  print("==============================")

  print(f"All Dice : {np.mean(all_d):.4f}")
  print(f"All IoU  : {np.mean(all_i):.4f}")

  print("------------------------------")

  print(f"Pos Dice : {np.mean(pos_d):.4f}")
  print(f"Pos IoU  : {np.mean(pos_i):.4f}")

  print("------------------------------")

  print(
    f"False Positive Rate : "
    f"{fp_count/(neg_count+1e-6):.4f}"
  )

  print("==============================")

  return np.mean(pos_d)

# =========================
# 5. Main（带负样本控制）
# =========================
def main():
  set_seed(42)

  if not torch.cuda.is_available():
    raise RuntimeError("CUDA required")

  device = torch.device("cuda:0")

  train_ds = BrainCTDataset("png/images/train", "png/masks/train")
  val_ds = BrainCTDataset("png/images/validation", "png/masks/validation")

  train_loader = DataLoader(train_ds, batch_size=8, shuffle=True)
  val_loader = DataLoader(val_ds, batch_size=1, shuffle=False)

  model = smp.Unet(
    encoder_name="resnet34",
    encoder_weights=None,
    in_channels=1,
    classes=1
  ).to(device)

  optimizer = optim.AdamW(model.parameters(), lr=1e-4)
  bce = nn.BCEWithLogitsLoss()
  focal = FocalTverskyLoss()
  dice_loss = smp.losses.DiceLoss(smp.losses.BINARY_MODE, from_logits=True)

  best_score = 0
  best_path = "best_model.pth"

  # =========================
  # Training
  # =========================
  for epoch in range(50):
    model.train()

    for img, mask, _ in train_loader:
      img, mask = img.to(device), mask.to(device)

      optimizer.zero_grad()
      out = model(img)

      loss = 0.5 * focal(out, mask) + 0.3 * bce(out, mask) + 0.2 * dice_loss(out, mask)
      loss.backward()
      optimizer.step()

    # validation
    score = evaluate(model, val_loader, device)

    if score > best_score:
      best_score = score
      torch.save(model.state_dict(), best_path)
      print(f"Saved best model, score={best_score:.4f}")

  # =========================
  # Final evaluation
  # =========================
  print("\nFinal evaluation with best model:")
  model.load_state_dict(torch.load(best_path, map_location=device))
  print("\nThreshold Search")
  for t in np.arange(0.40, 0.91, 0.05): # use 0.5
    print(f"\nThreshold = {t:.2f}")
    evaluate( model, val_loader, device, threshold=t, min_area=25 )

# =========================
if __name__ == "__main__":
  main()