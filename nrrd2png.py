import os
import glob
import numpy as np
import SimpleITK as sitk
import cv2
from tqdm import tqdm


def window_transform(img, wl=40, ww=110):
  """
  img: HU image
  return: uint8(0~255)
  """

  lower = wl - ww / 2
  upper = wl + ww / 2

  img = np.clip(img, lower, upper)

  img = (img - lower) / (upper - lower)

  img = (img * 255).astype(np.uint8)

  return img


# ==========================
# parameter
# ==========================
DIVISION = "validation"

INPUT_DIR = os.path.join(r"./NRRDs", DIVISION)

OUT_IMG_DIR = os.path.join(r"./png/images", DIVISION)
OUT_MASK_DIR = os.path.join(r"./png/masks", DIVISION)

os.makedirs(OUT_IMG_DIR, exist_ok=True)
os.makedirs(OUT_MASK_DIR, exist_ok=True)

WL = 40
WW = 110


# ==========================
# main
# ==========================
def main():
  orig_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*-orig.nrrd")))

  for orig_path in tqdm(orig_files):
    case_id = os.path.basename(orig_path).replace("-orig.nrrd", "")
    mask_path = os.path.join(INPUT_DIR, f"{case_id}-label_revise.nrrd")
     
    # read CT, shape: [z,h,w]
    ct_itk = sitk.ReadImage(orig_path)
    ct = sitk.GetArrayFromImage(ct_itk)
    
    # read Mask
    mask_itk = sitk.ReadImage(mask_path)
    mask = sitk.GetArrayFromImage(mask_itk)  
    assert ct.shape == mask.shape
    
    # for slice
    for z in range(mask.shape[0]):
      mask_slice = mask[z]    
      # if not np.any(mask_slice == 1): continue

      # CT -> window-level/window-width -> uint8
      ct_slice = ct[z]    
      ct_png = window_transform(ct_slice, wl=WL, ww=WW)
      
      mask_png = (mask_slice == 1).astype(np.uint8) # mask -> uint8

      img_name = f"{case_id}_{z:03d}.png"
      img_save_path = os.path.join(OUT_IMG_DIR, img_name)
      mask_save_path = os.path.join(OUT_MASK_DIR, img_name)
      cv2.imwrite(img_save_path, ct_png)
      cv2.imwrite(mask_save_path, mask_png)

  print("Done.")


# =========================
if __name__ == "__main__":
  main()