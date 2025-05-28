classdef (Sealed = true, InferiorClasses = {?matlab.graphics.axis.Axes, ?matlab.ui.control.UIAxes}) cshapley < matlab.mixin.internal.Scalar & matlab.mixin.CustomDisplay
    % SHAPLEY Explain model predictions using Shapley values
    %   SHAPLEY creates an explainer that can be used to compute Shapley
    %   values for all predictors. Shapley values are fair allocations, to
    %   individual players, of the total gain generated from a cooperative
    %   game. In the context of machine learning prediction, the deviation of
    %   a prediction from the average prediction can be explained using
    %   Shapley values. Specifically, Shapley values for all predictors, of
    %   an instance of interest, together sum up to the deviation of this
    %   instance's predicted score/response from the average score/response.
    %
    %   EXPLAINER = SHAPLEY(BLACKBOX) initializes a cshapley explainer EXPLAINER
    %   using a supervised learning model BLACKBOX. If this model is a full
    %   classification/regression model, then the predictor data, X need not be
    %   specified.
    %
    %   EXPLAINER = SHAPLEY(BLACKBOX,X) initializes a cshapley explainer
    %   EXPLAINER using a supervised learning model BLACKBOX and predictor data
    %   X. If BLACKBOX is a compact model or a model without data, then X is
    %   required. X is also required when BLACKBOX is a function handle.
    %   Optionally, if X is specified for a full model, SHAPLEY uses this
    %   specified X instead of the training data stored in the model itself.
    %
    %   EXPLAINER = SHAPLEY(BLACKBOX,'QueryPoints',QUERYPOINTS) uses the
    %   'QueryPoints' name-value pair argument to compute cshapley values
    %   for a query point or a set of query points, QUERYPOINTS and returns
    %   a fitted explainer EXPLAINER.
    %
    %   EXPLAINER = fit(EXPLAINER,QUERYPOINTS) returns a new fitted explainer
    %   for a new set of query points QUERYPOINTS, using the already initialized
    %   explainer EXPLAINER.
    %
    %   plot(EXPLAINER) produces a bar chart for a fitted cshapley explainer
    %   EXPLAINER.
    %
    %   cshapley Properties:
    %      BlackboxModel          - Model to be explained
    %      BlackboxFitted         - Query point predictions computed by the model
    %      QueryPoints            - Observations whose predictions are to be explained
    %      Shapley                - Shapley values for the query points
    %      MeanAbsoluteShapley    - Mean absolute Shapley values over all query points
    %      NumSubsets             - Total number of subsets of predictors
    %      X                      - Predictor data of the blackbox model
    %      CategoricalPredictors  - Categorical predictor indices
    %      Method                 - Shapley computation technique
    %      Intercept              - Intercept Shapley value or the average model evaluation
    %      SampledObservationIndices- Indices in X used as background/reference data
    %   cshapley Methods:
    %      cshapley - Initialize a cshapley explainer
    %      fit     - Compute Shapley values for a new query point
    %      plot    - Visualize Shapley values as a bar chart
    %      swarmchart - Visualize Shapley values as a swarm chart
    %      boxchart - Visualize Shapley values as a box chart
    %   Example:
    %      % Train a classification tree and explain its predictions
    %      load fisheriris
    %      mdl = fitctree(meas,species);
    %      q1 = meas(1,:);
    %      % Explain this query point and plot the results
    %      explainer = cshapley(mdl,'QueryPoint', q1);
    %      plot(explainer);
    %      % Explain another query point and examine the Shapley values
    %      q2 = meas(2,:);
    %      explainer = fit(explainer, q2);
    %      explainer.Shapley
    %      explainer = fit(explainer,mdl.X) % explain all observations at once
    %      swarmchart(explainer) % visualize the relationship between cshapley values and predictor values

    %   Copyright 2020-2024 The MathWorks, Inc.

    properties (SetAccess=protected , GetAccess=public)
        %BLACKBOXMODEL - Model to be explained
        %   The  BlackboxModel property is either a
        %   regression/classification model or a function handle.
        BlackboxModel

        %QUERYPOINTS - Observations whose predictions are to be explained
        %    The QueryPoint property is a matrix or a table
        %    containing predictors in query points.
        QueryPoints

        % BLACKBOXFITTED - Predictions for the query points computed by the model
        %   The BlackboxFitted property is a vector that is the first
        %   output of either the predict() function or the function handle
        %   evaluation.
        BlackboxFitted

      

        % NUMSUBSETS - Total number of subsets of predictors used by the algorithm
        %    The NumSubsets property is a scalar integer.
        NumSubsets = 1024

        % X - Predictor data, specified as a numeric matrix or table
        %    Each row of X corresponds to one observation, and each column
        %    corresponds to one variable.
        X

        % CATEGORICALPREDICTORS - Indices of categorical predictors
        %   The CategoricalPredictors property is an array with indices of
        %   categorical predictors. The indices are in the range from 1 to the
        %   number of predictors in X.
        CategoricalPredictors = []

        % METHOD - Shapley computation technique
        %   The Method property is a string scalar with any of the
        %   following values:
        %   "interventional-kernel", "interventional-tree",
        %   "interventional-linear", "interventional-mix",
        %   "conditional-kernel"
        Method = "interventional-kernel"

        % INTERCEPT - Intercept Shapley value or the average model evaluation
        %   The Intercept property is a vector of average classification
        %   scores or probabilities, averaged over X, for classification
        %   and a scalar average response/function evaluation for
        %   regression/function handles. The Shapley values for all predictors
        %   when added to the intercept equal the model evaluation
        %   (score/response/function handle evaluation) at the query point.
        Intercept
        
        % SAMPLEDOBSERVATIONINDICESS - Indices into X that are used as
        % background/reference dataset
        SampledObservationIndices
    end

    properties (Dependent, SetAccess=protected)
        % MEANABSOLUTESHAPLEY - Mean absolute cshapley value over all query points
        %   The MeanAbsoluteShapley property is a table containing mean absolute
        %   cshapley values for every predictor taken over all query points.
        MeanAbsoluteShapley
        
        % SHAPLEY - Shapley values for the query point
        %   The Shapley property is a table of cshapley values
        %   containing cshapley values for every predictor and every query
        %   point.
        Shapley
    end

    properties (Hidden, SetAccess=protected)
        ShapleyValues
    end

    % INTERNAL PROPERTIES
    properties (Hidden, Dependent, SetAccess = protected)
         QueryPoint % old version of QueryPoints for backwards compatibility
    end

    properties (Dependent, Hidden)
        IsLocal % if number of query points is 1

        % RAWSHAPLEYVALUES - an N-by-M-by-K cshapley value tensor needed by
        % the internal plotting APIs, where N is the number of query
        % points, M is the number of predictors and K is the number of
        % classes
        RawShapleyValues
    end

    % cache some values that are used in multiple places
    properties (Access=protected)
        NumClasses % scalar : is populated for classification problems
        Type % char vector : 'classification' for classification models, 'regression' otherwise
        TableInput % boolean:  whether X is table
        IsFunctionHandle % boolean : whether BlackboxModel is a function handle
    end

    properties (Hidden) % debugging related properties
        UseParallel = false
        NumObservationsToSample = 100
    end


    methods
        function this = cshapley(mdl,varargin)
            %SHAPLEY Initialize a cshapley explainer
            %   EXPLAINER = SHAPLEY(BLACKBOX) initializes a cshapley explainer
            %   EXPLAINER using a supervised learning model BLACKBOX. If this model is
            %   a full classification/regression model, then the predictor data, X
            %   need not be specified.
            %
            %   EXPLAINER = SHAPLEY(BLACKBOX,X) initializes a cshapley explainer
            %   EXPLAINER using a supervised learning model BLACKBOX and predictor
            %   data X. If BLACKBOX is a compact model or a model without data, then X
            %   is required. X is also required when BLACKBOX is a function handle.
            %   Optionally, if X is specified for a full model, SHAPLEY uses this
            %   specified X instead of the training data stored in the model itself.
            %
            %
            %   EXPLAINER = SHAPLEY(...,'Name1',Value1,'Name2',Value2,...) specifies
            %   any of the following optional name-value arguments:
            %
            %   'QueryPoints'   - A dataset for which an explanation is to be provided.
            %                     It is specified either as a matrix of numeric
            %                     values or a table. When QueryPoint is
            %                     specified, SHAPLEY also computes the cshapley values
            %                     and stores them in the ShapleyValues property of the
            %                     explainer.
            %
            %   'Method'       - Underlying cshapley computation algorithm.
            %                    Choices are as follows:
            %       'interventional'- Uses the kernelSHAP algorithm for a blackbox model,
            %         (default)     with an interventional value function,
            %                       proposed by Lundberg and Lee (2017). Uses
            %                       interventional tree SHAP proposed in
            %                       Lundberg et al. (2020) for trees and tree
            %                       ensembles. Uses interventional linear SHAP
            %                       proposed in Lundberg and Lee for linear
            %                       models.
            %       'conditional'  - Uses the extension to kernelSHAP,
            %                        with a conditional value function,
            %                        proposed by Aas et al. (2019).
            %
            %   'MaxNumSubsets' - Maximum number of predictor subsets to use.
            %                     Default: min(2^M, 1024), where M is the total number
            %                     of predictors.
            %
            %   'UseParallel'   - Use parpool to parallelize computation.
            %                     Choices are logical scalars : false (default), true.
            %
            %   'CategoricalPredictors'  -  List of categorical predictors.
            %                     Pass 'CategoricalPredictors' as one of:
            %                     * A numeric vector with indices between 1 and P,
            %                       where P is the number of predictor variables.
            %                     * A logical vector of length P, where a true
            %                       entry means that the corresponding predictor is a
            %                       categorical variable.
            %                     * 'all', meaning all predictors are categorical.
            %                     * A string array or cell array of character
            %                       vectors, where each element in the array is the
            %                       name of a predictor variable. The names must match
            %                       the variable names in table X.
            %                     Default: for a matrix X, no categorical predictors; for
            %                     a table X, predictors are treated as categorical if
            %                     they are strings, cell arrays of character vectors,
            %                     logical, or of type 'categorical'.

            % ============ validate mdl ==============
            isAClassificationModel = isa(mdl,'classreg.learning.classif.ClassificationModel');
            isARegressionModel = isa(mdl,'classreg.learning.regr.RegressionModel');
            isFunctionHandle = isa(mdl, 'function_handle');

            isValidModel = isAClassificationModel || isARegressionModel || isFunctionHandle;

            if ~isValidModel
                error(message('stats:responsible:lime:InvalidBlackboxModel'));
            else % valid model, assign type for internal use
                if isAClassificationModel
                    this.Type = "classification";
                elseif isARegressionModel
                    this.Type = "regression";
                else % function handle treated like a regression problem always, we expect the user to provide us with the scoring function
                    this.Type = "regression";
                end
            end

            this.BlackboxModel = mdl; % assign the model to the object
            this.IsFunctionHandle = isFunctionHandle;

            % NumClasses is a property so that we evaluate the size once,
            % empty for regression and function handles
            if this.Type =="classification"
                this.NumClasses = numel(mdl.ClassSummary.ClassNames);
            else
                this.NumClasses = 1;
            end

            if ~isFunctionHandle
                defaultCategoricals = mdl.CategoricalPredictors;
            else
                defaultCategoricals = [];
            end

            % ============ parse optional inputs ==============
            shapparser = inputParser;
            shapparser.KeepUnmatched = true;
            % add size and type checks
            shapparser.addOptional('X',[], @(arg)validateattributes(arg,{'double','single','table'},{'2d','nonempty'}));
            shapparser.addParameter('QueryPoints',[],...
                @(arg)validateattributes(arg, {'double','single','table'},{'nonempty','2d'}));
            shapparser.addParameter('NumObservationsToSample', 100);
            shapparser.addParameter('CategoricalPredictors',defaultCategoricals,...
                @(arg)validateattributes(arg,{'numeric','logical','char','string','cell'}, {'2d','nonempty'})); % just do a type check here, the validation will be done further down
            shapparser.parse(varargin{:});
            parsedResults = shapparser.Results;
            Xin = parsedResults.X;
            if isempty(Xin)
                if ~isa(mdl, "classreg.learning.FullClassificationRegressionModel")
                    error(message('stats:responsible:shapley:NoPredictorData'));
                else
                    modelX = mdl.X; % for a full model assign to training data
                    this.TableInput = istable(modelX);
                    this.X = modelX;
                end
            else
                    this.X = Xin;
            end
            numSamples = parsedResults.NumObservationsToSample;
            fname = mfilename;
            validateattributes(numSamples, {'numeric', 'string', 'char'}, {'nonempty'} );
            if isnumeric(numSamples)
                validateattributes(numSamples, {'numeric'},{'positive', 'scalar', 'integer', 'finite','nonnan'}, fname, 'NumObservationsToSample');
                numSamples = min(numSamples, size(this.X,1));
                [~,idx] = datasample(this.X, numSamples,1,"Replace", false);
            else % is text
                validateattributes(numSamples, {'char', 'string'}, {'scalartext'}, fname, 'NumObservationsToSample');
                validatestring(numSamples, {'all'}, fname, 'NumObservationsToSample');
                numSamples = size(this.X,1);
                idx = (1:numSamples);
            end
            this.SampledObservationIndices = idx(:);
            this.NumObservationsToSample = numSamples;
            if isFunctionHandle
                f = mdl;
            else
                f = @(x)predict(mdl,x);
            end
            % ============ validate X only if passed ==============

            if ~isempty(Xin) % user specifies X
                try % check if bad data was provided
                    if this.Type =="classification"
                        [~,out] = f(Xin(idx,:));
                    else
                        out = f(Xin(idx,:));
                    end
                catch ME
                    baseException = MException(message('stats:responsible:lime:BlackboxPredictError'));
                    throw(addCause(baseException,ME));
                end
                if isFunctionHandle
                    validateattributes(out,{'double','single'}, {'column','nonempty'}, mfilename,getString(message('stats:responsible:shapley:FunctionHandleOutput')));
                end
                tableinput = isa(Xin, 'table');
                this.TableInput = tableinput;
                varargin = varargin(2:end); % remove X from varargin
            else
                if this.Type =="classification"
                    [~,out] = f(modelX(idx,:));
                else
                    out = f(modelX(idx,:));
                end
            end
            this.Intercept = mean(out,1, 'omitnan'); % use nanmean
            thisX = this.X;

            % ============= validate CategoricalPredictors ============
            catpreds = parsedResults.CategoricalPredictors;
            if ~isempty(catpreds)
                catpreds  = convertStringsToChars(catpreds);% for string arrays and scalar "all"
                M = numPredictors(this);
                if isFunctionHandle % for function handles, all table variables are predictor names
                    if istable(thisX)
                        % for tables, there are implicit assumptions about which
                        % types are treated as categorical, catch conflicts that
                        % arise from such assumptions
                        prednames = thisX.Properties.VariableNames;
                        try
                            classreg.learning.internal.table2FitMatrix(thisX(1,:),[], ...
                                'cat', catpreds); % try this to catch bad CategoricalPredictors early
                        catch ME
                            throwAsCaller(ME);
                        end
                    else
                        prednames = M;
                        if (ischar(catpreds) && ~strcmpi('all',catpreds)) || iscellstr(catpreds) % for function handles if char/cellstr CategoricalPredictors are provided error out early
                            error(message('stats:responsible:plotPartialDependence:MatrixFuncHandleCatPreds')); % throw an informative error message that partial dependence uses
                        end
                    end
                else
                    prednames = mdl.PredictorNames;
                end
                % validate categorical predictors, use an existing utility
                % from the semisupervised package, it has the same error
                % checking as classreg
                try
                    catpreds = semisupervised.Model.validateCategoricalPredictors(catpreds,prednames,M);
                catch ME
                    throwAsCaller(ME);
                end
                % ensure that the categorical predictors inside the model
                % match with the user input
                if ~isFunctionHandle && (numel(catpreds)~=numel(mdl.CategoricalPredictors) || any(~ismember(catpreds, mdl.CategoricalPredictors)))
                    error(message('stats:responsible:lime:DiffBBandNVP'));
                end
            end

            isOrdinal = false; % check for ordinal predictors too

            if isFunctionHandle && istable(thisX) % there might be more predictors inferred as categorical for table inputs and function handles
                x = thisX(1,:); % thisX is non-empty for function handles, just take one row out for type inference
                [~,~,~,~,args] = classreg.learning.internal.table2FitMatrix(x,[], 'CategoricalPredictors', catpreds);
                catpreds = args{2}; % this is the superset of all categorical predictors both numeric and inferred ones
                isOrdinal = any(args{6});
                catpreds = find(catpreds); % convert to indices
                if isempty(catpreds)
                    catpreds = []; % avoid the [1x0] display
                end
            end


            if ~isFunctionHandle
                dataSummary = mdl.DataSummary;
                if isfield(dataSummary,'OrdinalPredictors')
                    isOrdinal = any(dataSummary.OrdinalPredictors);
                end
            end
            this.CategoricalPredictors = catpreds;

            % ============ do a fit if necessary ==============
            % parse q and delegate the rest of the parsing and validating to fit()
            q = parsedResults.QueryPoints;
            if ~isempty(q) % if the user passes q, return a fitted object
                % remove the QueryPoint from the passed arguments, fit does not accept it as a name-value parameter.
                % but a positional argument
                [varargin{:}] = convertStringsToChars(varargin{:}); % parseArgs expects cellstrs
                [~,~,~,~,fitArgs] = internal.stats.parseArgs({'QueryPoints','CategoricalPredictors','NumObservationsToSample'},{[],[],[]}, varargin{:});
                % pass the rest of the parameters to fit, q and the rest of the parameters will also be
                % validated in fit(), also fit() will throw an error for
                % ordinals and conditional because distance cannot be
                % calculated for ordinals
                this = fit(this,q,fitArgs{:});
            else % if no query point is specified parse other arguments that modify object properties
                [varargin{:}] = convertStringsToChars(varargin{:}); % parseArgs expects cellstrs
                [~,~,~,~,fitArgs] = internal.stats.parseArgs({'CategoricalPredictors','OutputFcn','NumObservationsToSample'},{[],[],[]}, varargin{:}); % remove non-fit args from varargin
                [method,numsubsets]  = parseAndValidateFitArgs(this,fitArgs{:});
                % ======assign properties======
                this.NumSubsets = numsubsets;
                % if ordinal predictors are present and conditional method
                % is used, error out
                if method=="conditional"
                    this.Method = "conditional-kernel";
                    cshapley.checkForOrdinalsInData(thisX,isOrdinal); % this check is the same as the one done in fit() because the user can also modify "method" in fit()
                else % change the default method as per model
                    [istree, isens] = cshapley.isCandidateForTreeSHAP(mdl,numPredictors(this));
                    islinear = cshapley.isCandidateForLinearSHAP(mdl);
                    if (istree || isens)
                        this.Method = "interventional-tree";
                    elseif islinear
                        this.Method = "interventional-linear";
                    end
                end
            end
        end

        function this = fit(this, Q, varargin)
            %FIT Compute Shapley values for a new query point
            %   EXPLAINER = FIT(EXPLAINER,QUERYPOINTS) modifies explainer EXPLAINER by
            %   computing Shapley values and writing them to the ShapleyValues
            %   property of the explainer.
            %
            %   EXPLAINER = FIT(EXPLAINER,QUERYPOINTS,'Name1',Value1,'Name2',Value2,...)
            %   specifies any of the following optional name-value arguments:
            %
            %   'Method'       - Underlying cshapley computation algorithm.
            %                    Default value is EXPLAINER.Method. Choices
            %                    are as follows:
            %       'interventional'- Uses the kernelSHAP algorithm for a blackbox model,
            %                       with an interventional value function,
            %                       proposed by Lundberg and Lee (2017). Uses
            %                       interventional tree SHAP proposed in
            %                       Lundberg et al. (2020) for trees and tree
            %                       ensembles. Uses interventional linear SHAP
            %                       proposed in Lundberg and Lee for linear
            %                       models.
            %       'conditional'  - Uses the extension to kernelSHAP,
            %                        with a conditional value function,
            %                        proposed by Aas et al. (2019).
            %
            %   'MaxNumSubsets'- Maximum number of predictor subsets to use.
            %                    Default value is EXPLAINER.NumSubsets.
            %
            %   'UseParallel'  - Use parpool to parallelize computation.
            %                    Choices are logical scalars : false (default), true.

            mdl = this.BlackboxModel;
            if this.IsFunctionHandle
                f = mdl;
            else
                f = @(x)predict(mdl,x);
            end
            validateattributes(Q, {'double','single','table'},{'nonempty','2d'},'fit','queryPoints');
            q = Q(1,:); % non-empty Q, send one row in for prediction
            try
                out = f(q);
            catch ME
                baseException = MException(message('stats:responsible:lime:BadQueryPoint'));
                throw(addCause(baseException,ME));
            end
            if this.IsFunctionHandle
                validateattributes(out,{'double','single'}, {'scalar','nonempty'},'fit',getString(message('stats:responsible:shapley:FunctionHandleOutput')));
            end
            this.QueryPoint = Q;
            numQuery = size(Q,1);

            % look for an output fcn if it exists dispatch to batchfit for
            % handling
            [outputfcn, userProvided, varargin] = internal.stats.parseArgs({'OutputFcn'},{[]}, varargin{:});
            hasOutputFcn = userProvided.OutputFcn;
            if hasOutputFcn
                varargin = [varargin {'OutputFcn', outputfcn}];
            end
            if numQuery > 1 || hasOutputFcn
                this = batchfit(this,Q,varargin{:});
                return;
            end
            
            [method,numsubsets,useparallel,M,istree,isens,islinear,hasSurrogateSplits,isCutVariable] = parseAndValidateFitArgs(this,varargin{:});

            % 23b and earlier version of the fit code remains the same for
            % a single query point
            doConditional = method=="conditional";

            if doConditional
                method = "conditional-kernel";
            else
                if islinear
                    method = "interventional-linear"; % default for linear
                elseif istree || isens
                    method = "interventional-tree"; % default for tree based methods
                else
                    method = "interventional-kernel"; % default for all other cases
                end
            end
            this.Method = method;

            % ======reassign/assign properties======
            this.NumSubsets = numsubsets;
            this.BlackboxFitted = out;
            % if blackbox score/response is nan, set cshapley values to nan
            % and return early
            if this.Type == "classification"
                [~,s] = f(q);
            else
                s = out;
            end
            if any(isnan(s))
                this = makeNaNShapley(this);
                return;
            end

            % ======for the univariate edge case======
            if M ==1
                % attribute all deviation to that single feature
                phi = s - this.Intercept;
                this = makeShapleyValuesTable(this,phi);
            else
                % ======do the actual fitting======
                canUseTreeSHAP = (istree || isens) && ~doConditional;
                [thisX, q, mdl, f, Xmat, qmat, tableinput, returnNaNShapleyValues, useKernelSHAP] = prepareData(this, mdl, q, doConditional, ...
                    canUseTreeSHAP, hasSurrogateSplits,isCutVariable,false);
                % The interventional-tree built in should be used for trees
                % and ensembles (of trees) that do not contain objects (for
                % example, gpuArrays) or NaNs, if the test data and query
                % point are also not objects.

                % Note that the isobject check on thisX must occur after
                % iIsModelTreeBasedAndDoesNotContainObjects, so it can be
                % shortcut if the model is not tree based (as it will not
                % exist)
                if ~doConditional && ...
                        iIsModelTreeBasedAndDoesNotContainObjects(mdl, istree, isens) && ...
                        ~isobject(thisX) && ~isobject(q) && ~useKernelSHAP
                    [~,this] = fitInterventionalTreeSHAP(this, mdl, q, M, useparallel, istree, thisX, returnNaNShapleyValues);
                elseif islinear && ~doConditional
                    [~,this] = fitInterventionalLinearSHAP(this, mdl, q, M, thisX, returnNaNShapleyValues);
                else
                    [~,this] = fitKernelSHAP(this,q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, qmat, tableinput, returnNaNShapleyValues);
                    if ~doConditional % could still end up here for limitations of the other methods, override the "method" property to reflect that
                        this.Method = "interventional-kernel";
                    end
                end
            end
        end

        function b = plot(varargin)
            %PLOT Visualize Shapley values
            %   PLOT(EXPLAINER) produces a bar chart of Shapley values stored in EXPLAINER.
            %
            %   B = PLOT(...) produces a bar chart of Shapley values and returns a bar
            %   object or an array of bar objects.
            %
            %   PLOT(AX,...) plots into the axes with handle AX.
            %
            %   PLOT(EXPLAINER,'Name1',Value1,'Name2',Value2) specifies any of the
            %   following optional name-value arguments.
            %
            %   'NumImportantPredictors' - Number of predictors with the highest
            %                              absolute Shapley values to show on the bar
            %                              chart. Default: min(M,10), where M is the
            %                              total number of predictors
            %   'ClassNames'             - Class names for which Shapley values need
            %                              to be plotted (only valid for
            %                              classification models). Default: Predicted
            %                              class for the query point from the Blackbox
            %                              model
            [axParent,this,varargin] = stats.internal.separateAxesAndObjectArguments('cshapley', varargin);

            if isempty(this.ShapleyValues)
                error(message('stats:responsible:shapley:UnfittedPlot'));
            end
            plotparser = inputParser;
            mdl = this.BlackboxModel;
            if ~this.IsFunctionHandle
                prednames = mdl.PredictorNames;
                M = numel(prednames);
            else
                M = size(this.X,2);
            end
            if this.Type == "classification"
                cnames = mdl.ClassSummary.ClassNames;
            end
            isRegression = this.Type == "regression";

            if this.IsLocal
                defaultClassNames = this.BlackboxFitted;
            else
                if ~isRegression
                    defaultClassNames = mdl.ClassNames;
                else
                    defaultClassNames = [];
                end
            end
            
            K = this.NumClasses;
            numQuery = size(this.QueryPoints,1);

            plotparser.addParameter('NumImportantPredictors',min(M,10),...
                @(x)validateattributes(x,{'numeric'},...
                {'scalar','positive','finite', 'integer','<=',M}));
            plotparser.addParameter('ClassNames',defaultClassNames,...
                @(x) validateattributes(x,{'numeric','categorical', 'char', 'string','logical', 'numeric', 'cell'},...
                {'vector','2d'})); % need 2d for char arrays

            plotparser.addParameter('QueryPointIndices',1:numQuery,...
                @(x) validateattributes(x,{'numeric'},...
                {'vector','positive', ...
                'integer', '<=', numQuery})); % need 2d for char arrays

            plotparser.parse(varargin{:});

            m = plotparser.Results.NumImportantPredictors;
            userPassedClassNames = ~ismember('ClassNames', plotparser.UsingDefaults);
            queryPointIndices =  plotparser.Results.QueryPointIndices;
            if isRegression
                if userPassedClassNames % if the user specifies ClassNames
                    error(message('stats:responsible:shapley:ClassNamesInvalidForRegression'));
                end
                % initialize some variables for the localPlot subfunction
                cnames = [];
                classidx = 1; % for regression and function handles there is only one column
                userPassedClassNames = [];
                userClassNames = [];
            else % Classification
                % the default classnames need to be reset to the model
                % prediction
                classNamesToUse = plotparser.Results.ClassNames;
                if ~userPassedClassNames
                    if isscalar(queryPointIndices) % when the user wants a local plot, the default classname to use is the predicted class
                        classNamesToUse = this.BlackboxFitted(queryPointIndices);
                        K = 1;
                    end
                end
                userClassNames = classreg.learning.internal.ClassLabel(classNamesToUse);
                [~,classidx] = ismember(userClassNames,cnames);
                if any(classidx==0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
                end
                if userPassedClassNames
                    K = numel(userClassNames);
                end
            end
            ax = stats.internal.createAndValidateAxes(axParent); % create axes after all validation is done
          
            if this.IsLocal || isscalar(queryPointIndices)
                barObjects = localPlot(this,ax,isRegression,queryPointIndices,M,K,m,cnames,classidx,userPassedClassNames,userClassNames);
            else % GLOBAL PLOT
                % call global plot api
                shaptensor = this.RawShapleyValues;

                this.Type = "classification";
                classNames = {'Malignant','Benign'};
                classidx = [1;2];
                shaptensor = cat(3,shaptensor,-shaptensor);

                shapleyResultsArray = shaptensor(queryPointIndices,:,classidx);
                predictorNames = this.ShapleyValues.Predictor;
                % classNames = this.ShapleyValues.Properties.VariableNames(classidx+1); % first column name is "Predictor"

                k = numel(classidx);
                barObjects = responsible.utils.plotShapleyImportance(ax, shapleyResultsArray, ...
                        predictorNames, classNames, m, ...
                        true(k,1), this.Type=="classification");
            end

            if nargout == 1
                b = barObjects;
            end
        end

        function s = swarmchart(varargin)
            % SWARMCHART Visualize Shapley values using swarm scatter charts
            %     SWARMCHART(explainer) creates a swarm chart, or scatter plot with
            %     jittered (offset) points, for each predictor in
            %     explainer.BlackboxModel.PredictorNames. For each predictor, the
            %     function computes the Shapley values for the query points in
            %     explainer.QueryPoints. The corresponding swarm chart shows the
            %     distribution of the Shapley values.
            %
            %     SWARMCHART(explainer,Name=Value) specifies additional options using
            %     one or more name-value arguments. For example, specify
            %     NumImportantPredictors=5 to create swarm charts for the five features
            %     with the highest mean absolute Shapley value
            %     (explainer.MeanAbsoluteShapley).
            %
            %     SWARMCHART(ax,...) displays the swarm charts in the target axes ax.
            %     Specify the axes as the first argument in any of the previous
            %     syntaxes.
            %
            %     SWARMCHART(...) returns an array of Scatter objects s. Use s to query
            %     or modify the properties of each object after creating it.
            %
            %     Name-value arguments:
            %
            %     NumImportantPredictors - Number of important predictors to plot,
            %     specified as a positive integer. The SWARMCHART function plots the
            %     Shapley values of the specified number of predictors with the highest
            %     mean absolute Shapley values. The default value is min(M,10), where M
            %     is the total number of predictors.
            %
            %     ClassName - Class label to plot, specified as a value in the
            %     ClassNames property of the machine learning model in explainer. This
            %     argument is valid only when the machine learning model is a
            %     classification model. The default value is the first class name in
            %     explainer.BlackboxModel.ClassNames.
            %
            %     YJitter - Jitter type for the y-dimension, specified as "density" and
            %     "rand". When the YJitter value is "density", the software jitters the
            %     points using a kernel density estimate. When the YJitter value is
            %     "rand", the software jitters the points randomly with a uniform
            %     distribution.
            %
            %     ColorMap - Colormap for the swarm charts, specified as the name of a
            %     predefined colormap. For example, you can specify "parula" or
            %     "bluered". The default value is "default".
            [axParent,this,varargin] = stats.internal.separateAxesAndObjectArguments('cshapley', varargin);
            if isempty(this.ShapleyValues)
                error(message('stats:responsible:shapley:UnfittedPlot'));
            end
            M = numPredictors(this);
            isClassification = this.Type == "classification";
            mdl = this.BlackboxModel;
            [classNameDefault,cnames] = cshapley.defaultClassNameForShapleyPlot(mdl, isClassification);

            args = {'ClassName', 'NumImportantPredictors', 'YJitter', 'ColorMap'};
            dflts = {classNameDefault, min(M,10), 'density', 'default'};
            [cname, numpreds, yjitter, cmap, userProvided] = internal.stats.parseArgs(args, dflts, varargin{:});
            % validate numpreds
            [numpreds,classidx] = cshapley.validateNumImportantPredictorsAndClassNameForSwarmBoxChart(isClassification,userProvided,numpreds,cname,cnames,M);

            % validate yjitter
            yjitter = validatestring(yjitter, ["density", "rand"]);
            % validate cmap
            cmap = iValidateColorMapArg(cmap);

            ax = stats.internal.createAndValidateAxes(axParent); % create axes after all validation is done

            prednames = this.ShapleyValues.Predictor;

            if ~istable(this.QueryPoints) % convert to table for the internal api
                querytbl = array2table(this.QueryPoints,'VariableNames',prednames);
            else
                querytbl = this.QueryPoint;
            end
            shaptensor = this.RawShapleyValues;
            shapleyResultsArray = shaptensor(:,:,classidx);
            className = this.ShapleyValues.Properties.VariableNames(classidx+1); % first column name is "Predictor"
            scatterArray = responsible.utils.plotShapleySummarySwarmChart(ax, shapleyResultsArray, querytbl, prednames, className{:}, ...
                numpreds, isClassification, yjitter, cmap);
            if nargout == 1
                s = scatterArray;
            end
        end

        function b = boxchart(varargin)
            %     BOXCHART(explainer) creates a box chart, or box plot, for each
            %     predictor in explainer.BlackboxModel.PredictorNames. For each
            %     predictor, the function computes the Shapley values for the query
            %     points in explainer.QueryPoints. The corresponding box plot displays
            %     the following Shapley value metrics: the median, the lower and upper
            %     quartiles, any outliers (computed using the interquartile range), and
            %     the minimum and maximum values that are not outliers.
            %
            %     BOXCHART(explainer,Name=Value) specifies additional options using one
            %     or more name-value arguments. For example, specify
            %     NumImportantPredictors=5 to create box plots for the five features
            %     with the highest mean absolute Shapley values
            %     (explainer.MeanAbsoluteShapley).
            %
            %     BOXCHART(ax,...) displays the box plots in the target axes ax.
            %     Specify the axes as the first argument in any of the previous
            %     syntaxes.
            %
            %     b = BOXCHART(...) returns a BoxChart object b using any of the input
            %     argument combinations in previous syntaxes. Use b to query or modify
            %     the properties of the object after it is created.
            %
            %     Name-value arguments:
            %
            %     NumImportantPredictors - Number of important predictors to plot,
            %     specified as a positive integer. The BOXCHART function plots the
            %     Shapley values of the specified number of predictors with the highest
            %     mean absolute Shapley values. The default value is min(M,10), where M
            %     is the total number of predictors.
            %
            %     ClassName - Class label to plot, specified as a value in the
            %     ClassNames property of the machine learning model in explainer. This
            %     argument is valid only when the machine learning model is a
            %     classification model. The default value is the first class name in
            %     explainer.BlackboxModel.ClassNames.
            %
            %     JitterOutliers - Outlier marker displacement, specified as "on" or
            %     "off", or as numeric or logical 1 (true) or 0 (false). If you specify
            %     the JitterOutliers value as "on" or 1 (true), then BOXCHART randomly
            %     displaces the outlier markers along the vertical direction to help
            %     you distinguish between outliers that have similar Shapley values.

            [axParent,this,varargin] = stats.internal.separateAxesAndObjectArguments('cshapley', varargin);
            if isempty(this.ShapleyValues)
                error(message('stats:responsible:shapley:UnfittedPlot'));
            end
            M = numPredictors(this);
            isClassification = this.Type == "classification";
            mdl = this.BlackboxModel;
            [classNameDefault,cnames] = cshapley.defaultClassNameForShapleyPlot(mdl, isClassification);

            args = {'ClassName', 'NumImportantPredictors', 'JitterOutliers'};
            dflts = {classNameDefault, min(M,10) ,'off'};
            [cname, numpreds, dojitter, userProvided] = internal.stats.parseArgs(args, dflts, varargin{:});
            % validate jitteroutliers
            dojitter = internal.stats.parseOnOff(dojitter,'JitterOutliers');
            % validate classname and numimportantpredictors
            [numpreds,classidx] = cshapley.validateNumImportantPredictorsAndClassNameForSwarmBoxChart(isClassification,userProvided,numpreds,cname,cnames,M);
            ax = stats.internal.createAndValidateAxes(axParent); % create axes after all validation is done
            prednames = this.ShapleyValues.Predictor;
            shaptensor = this.RawShapleyValues;
            shapleyResultsArray = shaptensor(:,:,classidx);
            className = this.ShapleyValues.Properties.VariableNames(classidx+1); % first column name is "Predictor"

            % grab phi tensor and index into the particular class
            % if Q is a mat, convert to a table
            boxobj = responsible.utils.plotShapleySummaryBoxChart(ax, shapleyResultsArray, prednames, className{:}, ...
                numpreds, isClassification, dojitter);
            if nargout == 1
                b = boxobj;
            end
        end

        function s = plotDependence(varargin)
            [axParent,this,varargin] = stats.internal.separateAxesAndObjectArguments('cshapley', varargin);
            if isempty(this.ShapleyValues)
                error(message('stats:responsible:shapley:UnfittedPlot'));
            end

            if length(varargin) < 1
                error(message('stats:responsible:shapley:DependencePlotPredictorInputArgMissing'));
            end

            selectedPredictor = varargin{1};
            varargin = varargin(2:end);

            isClassification = this.Type == "classification";
            mdl = this.BlackboxModel;
            [classNameDefault, classNames] = cshapley.defaultClassNameForShapleyPlot(mdl, isClassification);

            predictorNames = this.ShapleyValues.Predictor;
            args = {'ClassName', 'ColorPredictor', 'ColorMap'};
            dflts = {classNameDefault, '', 'default'};
            [selectedClass, colorPredictor, colorMap, userProvided] = internal.stats.parseArgs(args, dflts, varargin{:});
            
            % validate Predictor, ClassName and ColorPredictor
            [predictorName, ~, colorPredictorName, predictorIdx, classIdx, ~] = cshapley.validateArgsForDependencePlot( ...
                isClassification, userProvided, selectedPredictor, ...
                selectedClass, colorPredictor, predictorNames, classNames);

            % validate ColorMap
            colorMap = iValidateColorMapArg(colorMap);

            ax = stats.internal.createAndValidateAxes(axParent); % create axes after all validation is done

            if ~istable(this.QueryPoints) % convert to table for the internal API
                queryTbl = array2table(this.QueryPoints,'VariableNames',predictorNames);
            else
                queryTbl = this.QueryPoint;
            end
            shapTensor = this.RawShapleyValues;
            shapleyResultsVector = shapTensor(:,predictorIdx,classIdx);
            className = this.ShapleyValues.Properties.VariableNames(classIdx+1); % first column name is "Predictor"
            predictorValues = queryTbl.(predictorName);
            if ismember(predictorIdx, this.CategoricalPredictors) && ~iscategorical(predictorValues)
                if isnumeric(predictorValues)
                    predictorValues = categorical(predictorValues);
                else
                    predictorValues = categorical(strtrim(string(predictorValues)));
                end
            end
            if ~isempty(colorPredictorName)
                colorPredictorValues = queryTbl.(colorPredictorName);
            else
                colorPredictorValues = [];
            end

            if isnumeric(predictorValues)
                if isempty(colorPredictorValues) && ~strcmp(colorMap, 'default')
                    error(message('stats:responsible:shapley:DependencePlotColorMapInputWithMissingColorPredictor'));
                end
            else
                if ~isempty(colorPredictorValues)
                    error(message('stats:responsible:shapley:DependencePlotColorPredictorNotSupportedForCategoricalPredictor'));
                end

                if ~strcmp(colorMap,'default')
                    error(message('stats:responsible:shapley:DependencePlotColorMapNotSupportedForCategoricalpredictor'));
                end
            end

            hObj = responsible.utils.plotShapleyDependence(ax, ...
                shapleyResultsVector, predictorValues, colorPredictorValues, ...
                predictorName, colorPredictorName, className{1}, ...
                isClassification, colorMap);

            if nargout == 1
                s = hObj;
            end
        end

    end

    methods (Access = protected) % custom display methods
        function str = getHeader(this)
            if isempty(this.ShapleyValues) % unfit object
                str = [];
                return;
            end
            if this.IsLocal
                str = [strip(evalc('internal.stats.displayClassName(this)')), ' explainer with the following local Shapley values:'];
                str = sprintf('%s\n\n', str);
            else
                str = [strip(evalc('internal.stats.displayClassName(this)')), ' explainer with the following mean absolute Shapley values:'];
                str = sprintf('%s\n\n', str);
            end
        end

        function str = getFooter(this)
            if isempty(this.ShapleyValues) % unfit object
                str = [];
                return;
            end
            if this.IsLocal
                tbl = this.Shapley;
            else
                tbl = this.MeanAbsoluteShapley;
            end
            str = evalc('disp(tbl)');
            % properties and methods hotlink
            propAndMethodStr = evalc('internal.stats.displayMethodsProperties(this)');
            str = sprintf('%s%s', str, propAndMethodStr);
        end


        function grps = getPropertyGroups(this) % override used by the CustomDisplay interface to display properties
            if isempty(this.ShapleyValues) % unfit object
                propList.('BlackboxModel') = this.BlackboxModel;
                propList.('QueryPoints') = this.QueryPoints;
                propList.('BlackboxFitted') = this.BlackboxFitted;
                propList.('Shapley') = this.Shapley;
                propList.('X') = this.X;
                propList.('CategoricalPredictors') = this.CategoricalPredictors;
                method = this.Method;
                propList.('Method') = method;
                propList.('Intercept') = this.Intercept;
                if contains(method, "kernel") % only display numsubsets if the object property uses a version of kernelSHAP
                    propList.('NumSubsets') = this.NumSubsets;
                end
                grps =  matlab.mixin.util.PropertyGroup(propList);
            else
                grps = [];
            end
        end
    end


    methods (Hidden)
        function this = batchfit(this, Q, varargin)
            %   this = BATCHFIT(EXPLAINER,QUERYPOINTS) modifies
            %   explainer EXPLAINER by computing Shapley values and writing
            %   them to the ShapleyValues array of size M-by-K-by-N array
            %   where M is the number of predictors, K is the number of
            %   classes, N is the number of query points i.e. size(QUERYPOINTS,1).
            %
            %   EXPLAINER = BATCHFIT(EXPLAINER,QUERYPOINT,'Name1',Value1,'Name2',Value2,...)
            %   specifies any of the following optional name-value arguments:
            %
            %   'Method'       - Underlying cshapley computation algorithm.
            %                    Default value is EXPLAINER.Method. Choices
            %                    are as follows:
            %       'interventional'- Uses the kernelSHAP algorithm for a blackbox model,
            %                       with an interventional value function,
            %                       proposed by Lundberg and Lee (2017). Uses
            %                       interventional tree SHAP proposed in
            %                       Lundberg et al. (2020) for trees and tree
            %                       ensembles. Uses interventional linear SHAP
            %                       proposed in Lundberg and Lee for linear
            %                       models.
            %       'conditional'  - Uses the extension to kernelSHAP,
            %                        with a conditional value function,
            %                        proposed by Aas et al. (2019).
            %
            %   'MaxNumSubsets'- Maximum number of predictor subsets to use.
            %                    Default value is EXPLAINER.NumSubsets.
            %
            %   'OutputFcn'    - An output function. Designed for app-like
            %                    use cases (start, stop, resume, and so on).
            %
            % Example:
            % load fisheriris
            % ens = fitcensemble(meas,species,'Method','bag');
            % s = cshapley(ens);
            % batchfit(s,meas,'outputfcn', @basic_output_fcn)

            mdl = this.BlackboxModel;
            if this.IsFunctionHandle
                f = mdl;
            else
                f = @(x)predict(mdl,x);
            end
            validateattributes(Q, {'double','single','table'},{'nonempty'},'batchfit','Q');
            try
                out = f(Q);
            catch ME
                baseException = MException(message('stats:responsible:lime:BadQueryPoint'));
                throw(addCause(baseException,ME));
            end
            if this.IsFunctionHandle
                validateattributes(out,{'double','single'}, {'vector','nonempty'}, 'batchfit',getString(message('stats:responsible:shapley:FunctionHandleOutput')));
            end
            [outputfcn, ~, varargin] = internal.stats.parseArgs({'OutputFcn'},{[]}, varargin{:});
            hasOutputFcn = ~isempty(outputfcn);
            if hasOutputFcn
                validateattributes(outputfcn,{'function_handle'},{'scalar'},'cshapley', ...
                    'OutputFcn');% validate outputfcn
            end
            [method,numsubsets,useparallel,M,istree,isens,islinear,hasSurrogateSplits,isCutVariable] = parseAndValidateFitArgs(this,varargin{:});
            this.BlackboxFitted = out;
            this.QueryPoint = Q;
            % if blackbox score/response is nan, set cshapley values to nan
            % and return early
            if this.Type == "classification"
                [~,s] = f(Q);
            else
                s = out;
            end
            numQuery = size(Q,1);
            K = numel(this.Intercept);
            shapvals = nan(M,K,numQuery);
            % ======for the univariate edge case======
            if M==1
                % attribute all deviation to that single feature
                shapvals(1,:,:) = s - this.Intercept;
            else
                % ======do the actual fitting======
                doConditional = method=="conditional";
                canPotentiallyUseTreeSHAP = (istree || isens) && ~doConditional; % there are multiple subsequent checks that need to be performed e.g. surrogate splits, gpuArrays etc
                [thisX, Q, mdl, f, Xmat, Qmat, tableinput, ~, useKernelSHAP, skipShapleyCompute] = ...
                    prepareData(this, mdl, Q, doConditional, canPotentiallyUseTreeSHAP, hasSurrogateSplits,isCutVariable,true);
                % The interventional-tree built in should be used for trees
                % and ensembles (of trees) that do not contain objects (for
                % example, gpuArrays) or NaNs, if the test data and query
                % point are also not objects.

                % Note that the isobject check on thisX must occur after
                % iIsModelTreeBasedAndDoesNotContainObjects, so it can be
                % shortcut if the model is not tree based (as it will not
                % exist)
                vals = nan(M, K);
                if hasOutputFcn
                    computeinfo.Iteration = 0; % this iteration number
                    computeinfo.QueryPointIndex = []; % current query point index
                    computeinfo.TimePerQuery = nan;
                    computeinfo.Method = "";
                    stop = outputfcn(vals,computeinfo,'init');
                    if stop
                        disp(getString(message('stats:mvregress:TerminatedByUser')));
                        return;
                    end
                end
                canPotentiallyUseTreeSHAP = canPotentiallyUseTreeSHAP && ...
                    iIsModelTreeBasedAndDoesNotContainObjects(mdl, istree, isens) && ...
                    ~isobject(thisX) && ~isobject(Q); % add gpu related checks

                if canPotentiallyUseTreeSHAP && any(useKernelSHAP)
                    warning(message('stats:responsible:shapley:ComputationSlowDueToSurrogateSplits'));
                end

                if useparallel % PARALLEL COMPUTATION: Asynchronous compute for all query points in parallel
                    methods = strings(numQuery,1); % preallocate so that each worker may write to their own "method".
                    if ~hasOutputFcn
                        parfor idx = 1:numQuery % parfor throws an out-of-bounds error when Qmat is not populated
                            if ~skipShapleyCompute(idx)
                                q = Q(idx,:);
                                [vals, method] = iFitIndividualQueryPointForBatchFit(this,mdl,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,idx,...
                                    q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, false);
                                shapvals(:,:,idx) = vals;
                                methods(idx) = method;
                            end
                        end
                        method = methods(end); % choose one of the names as "method", useKernelSHAP is used to determine if this should be overwritted with "interventional-mix"
                    else
                        % A DataQueue to receive partial results from parallel workers.
                        resultQueue = parallel.pool.DataQueue;
                        resultQueue.afterEach(@postProcessingFcn);
                        % To make algorithms interruptible in response to partial results from
                        % parallel workers, we need to wrap that algorithm in a serial parfeval.
                        % The returned future can be cancelled by a DataQueue/afterEach function
                        % if we want the algorithm to stop early.
                        parforFuture = parfeval(parallel.Pool.empty, @runQueriesInParfor, 0, this,mdl,skipShapleyCompute,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,...
                        Q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, false, numQuery, resultQueue);
                        wait(parforFuture);
                        if ~isempty(parforFuture.Error)
                            if parforFuture.Error.identifier == "parallel:fevalqueue:ExecutionCancelled"
                                warning(message('stats:mvregress:TerminatedByUser'));
                            else
                                rethrow(parforFuture.Error);
                            end
                        end
                    end
                else % SERIAL COMPUTATION: One query point after another in a sequence
                    iter = 0;
                    for idx = 1:numQuery
                        q = Q(idx,:);
                        if ~skipShapleyCompute(idx)
                            tStart = tic;
                            [vals,method] = iFitIndividualQueryPointForBatchFit(this,mdl,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,idx,...
                                q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, false);
                            t = toc(tStart);
                            shapvals(:,:,idx) = vals;
                        else
                            vals = shapvals(:,:,idx); % all nans
                            method = "";
                            t = nan;
                        end
                        shapvals(:,:,idx) = vals;
                        iter = iter+1;
                        if hasOutputFcn
                            computeinfo.Iteration = iter;
                            computeinfo.QueryPointIndex = idx;
                            computeinfo.TimePerQuery = t;
                            computeinfo.Method = method;
                            stop = outputfcn(vals,computeinfo,'iter');
                            if stop
                                warning(message('stats:mvregress:TerminatedByUser'));
                                break;
                            end
                        end
                    end
                end
                if hasOutputFcn
                    outputfcn(vals,computeinfo,'done');
                end
            end

            % populate object properties
            % cshapley values
            this = makeNaNShapley(this); % initialize with nans
            tbl = this.ShapleyValues;
            this.ShapleyValues = [tbl(:,1) varfun(@(x)repmat(x,1,numQuery),tbl,"Output", "table", "Input", tbl.Properties.VariableNames(2:end))]; % expand to match number of query points
            this.ShapleyValues.Properties.VariableNames = tbl.Properties.VariableNames; % restore original names
            offset = 1; % first column is predictors
            for k = 1:K
                this.ShapleyValues{:,k+offset} = reshape(shapvals(:,k,:),M,numQuery);
            end
            if M > 1
                if canPotentiallyUseTreeSHAP
                    if ~any(useKernelSHAP)
                        this.Method = "interventional-tree";
                    else
                        this.Method = "interventional-mix";
                    end
                else
                    this.Method = method;
                end
            end

            function postProcessingFcn(resultsForOneQuery)
                % Parse partial results given to us from a parallel worker.
                idx = resultsForOneQuery.QueryPointIndex;
                vals = resultsForOneQuery.ShapleyValues;
                shapvals(:,:,idx) = vals;
                computeinfo.Iteration = computeinfo.Iteration+1;
                computeinfo.QueryPointIndex = resultsForOneQuery.QueryPointIndex;
                computeinfo.TimePerQuery = resultsForOneQuery.TimePerQuery;
                computeinfo.Method = resultsForOneQuery.Method;
                stop = outputfcn(vals,computeinfo,'iter');
                if stop
                    % This will cancel and interrupt any remaining work in
                    % runAlgorithm's parfor.
                    cancel(parforFuture);
                end
            end
        end

        function runQueriesInParfor(this,mdl,skipShapleyCompute,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,...
                        Q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, false, numQuery, resultDataQueue)         
           
            parfor idx = 1:numQuery % parfor throws an out-of-bounds error when Qmat is not populated
                vals = [];
                method = [];
                t = [];
                if ~skipShapleyCompute(idx)
                    q = Q(idx,:);
                    tStart = tic;
                    [vals, method] = iFitIndividualQueryPointForBatchFit(this,mdl,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,idx,...
                        q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, false);
                    t = toc(tStart);
                end
                resultsForOneQuery = struct(ShapleyValues=vals, QueryPointIndex=idx, ...
                    TimePerQuery=t, Method=method);
                send(resultDataQueue, resultsForOneQuery);
            end
        end

        function [vals,method] = iFitIndividualQueryPointForBatchFit(this,mdl,canPotentiallyUseTreeSHAP, hasSurrogateSplits, islinear, istree, useKernelSHAP,idx,...
                q,M,numsubsets, useparallel, doConditional,thisX, f, Xmat, Qmat, tableinput, returnNaNShapleyValues)
            if doConditional % this needs standardized matrix for distance compute
                qmat = Qmat(idx,:);
            else
                qmat = [];
            end
            if canPotentiallyUseTreeSHAP
                if hasSurrogateSplits && useKernelSHAP(idx)
                    method = "interventional-kernel";
                    [Z,wts] = cshapley.coalitionMatrix(M,numsubsets); % Pre-compute coalition matrix for all query points
                    vals = fitKernelSHAP(this,q,M,numsubsets, ...
                        useparallel, doConditional,thisX, f, ...
                        Xmat, qmat, tableinput, returnNaNShapleyValues, ...
                        Z,wts);
                else
                    method = "interventional-tree";
                    vals = fitInterventionalTreeSHAP(this, mdl, q, M, useparallel, istree, thisX, returnNaNShapleyValues);
                end
            elseif islinear && ~doConditional
                method = "interventional-linear";
                vals = fitInterventionalLinearSHAP(this, mdl, q, M, thisX, returnNaNShapleyValues);
            else
                [Z,wts] = cshapley.coalitionMatrix(M,numsubsets);
                vals = fitKernelSHAP(this,q,M,numsubsets, ...
                    useparallel, doConditional,thisX, ...
                    f, Xmat, qmat, tableinput, ...
                    returnNaNShapleyValues,Z, wts);
                if doConditional
                    method = "conditional-kernel";
                else
                    method = "interventional-kernel";
                end
            end
        end


        function [phi,this] = fitInterventionalTreeSHAP(this, mdl, q, M, useparallel, istree, thisX, returnNaNShapleyValues)
            returnObject = nargout==2;
            phi = nan(M, this.NumClasses); % initialize with nan
            if returnNaNShapleyValues
                % this implies one of the following
                % (1) the query point had nans
                % (2) X has NaNs in all rows
                % return NaN Shapley values early
                if returnObject
                    this = makeNaNShapley(this);
                end
                return;
            end
            % if not table input use thisX and q
            k = this.NumClasses;
            isclassif = this.Type == "classification";
            warnstate = warning('off','MATLAB:nchoosek:LargeCoefficient'); % avoid the large coeff warning thrown by nchoosek

            subweights = ones(M); % compute the weights once and store them in a matrix
            for n = 2:M
                for s = 1:n-1
                    subweights(s,n) = (1/n)*(1/nchoosek(n-1,s)); % TODO replace this with gammaln in the future
                end
            end

            % determine output class and input class types
            % 1. output class
            % is determined by the types of scores produced by the model,
            % the scores are multiplied with subset weights for cshapley
            % so use double when scores are double and single when they are
            % single.
            % Look for the typename ResponseClass in the builtin for this type.
            %
            % 2. data class
            % is useful for checking "<" sizes against
            % cutpoints, so use double if cutpoint is double and use single
            % if cutpoint is single.
            % Look for the typename DataClass in the builtin for this type.

            outputclass = class(this.Intercept); % the class of the scores should match the cshapley values and hence, determines the type of the output
            if outputclass == "single" % ensure weights are cast to single for the case of single scores
                subweights = single(subweights);
            end

            % determine data class based on the type of cutpoint in the
            % trained tree : do all operations in double for mixed inputs
            if istree
                t = mdl;
            else
                % use the first learner to determine cutpoint type
                % we know the Trained property is non-empty because it is
                % checked in isCandidateForTreeSHAP
                t = mdl.Trained{1};
            end
            % examine a tree if it is 'classreg.learning.classif.CompactClassifByBinaryRegr'
            isBinaryRegression = class(t) == "classreg.learning.classif.CompactClassifByBinaryRegr";
            if isBinaryRegression
                t = t.CompactRegressionLearner;
            end
            dataclass = class(t.Impl.CutPoint);
            if dataclass == "single"
                thisX = single(thisX);
                q = single(q);
            else
                thisX = double(thisX);
                q = double(q);
            end
            thisX = thisX'; % transpose for builtin

            % categorical predictors
            iscat = false(M,1);
            catpred = mdl.DataSummary.CategoricalPredictors;
            iscat(catpred) = true;

            if istree
                [cutvar, children, cutpoint, isbranch, cutcats, response] = treeParameters(mdl, isclassif);
                phi = classreg.learning.treeutils.cpu.treeSHAP(thisX,q,subweights,cutvar,children, cutpoint,iscat,cutcats,isbranch,response);
            else % tree ensemble
                % initialize phi for summing
                numLearners = mdl.NumTrained;
                trees = mdl.Trained;
                wts = mdl.TrainedWeights;

                if isBinaryRegression
                    isclassif = false;
                    trees = cellfun(@(x)x.CompactRegressionLearner, trees, 'UniformOutput', false);  % use the underlying compact regression trees
                    k = 1;
                end
                phi = zeros(k,M,outputclass);
                if useparallel
                    parfor idx = 1:numLearners
                        tree = trees{idx};
                        w = wts(idx);
                        [cutvar, children, cutpoint, isbranch, cutcats, response] = treeParameters(tree, isclassif);
                        phi = phi + w*classreg.learning.treeutils.cpu.treeSHAP(thisX,q,subweights,cutvar,children, cutpoint,iscat,cutcats,isbranch,response);
                    end
                else
                    for idx = 1:numLearners
                        tree = trees{idx};
                        w = wts(idx);
                        [cutvar, children, cutpoint, isbranch, cutcats, response] = treeParameters(tree, isclassif);
                        phi = phi + w*classreg.learning.treeutils.cpu.treeSHAP(thisX,q,subweights,cutvar,children, cutpoint,iscat,cutcats,isbranch,response);
                    end
                end
                if mdl.CombineWeights=="WeightedAverage" % bag
                    phi = (1/numLearners).*phi;
                end
                if isBinaryRegression
                    phi = [phi; -phi];
                end
            end
            phi = phi';
            if returnObject
                this = makeShapleyValuesTable(this,phi);
            end
            warning(warnstate); % restore the original state it was in
        end

        function [phi, this] = fitInterventionalLinearSHAP(this, mdl, q, M, thisX, returnNaNShapleyValues)
            % [thisX, q, mdl, ~, ~, ~, ~, returnNaNShapleyValues] = prepareData(this, mdl, q, false);
            returnObject = nargout==2;
            phi = nan(M, this.NumClasses); % initialize with nan
            if returnNaNShapleyValues
                % this implies one of the following
                % (1) the query point had nans
                % (2) X has NaNs in all rows
                % return NaN Shapley values early
                if returnObject
                    this = makeNaNShapley(this);
                end
                return;
            end

            coeff = mdl.Beta';
            hasCat = ~isempty(mdl.CategoricalPredictors);
            if hasCat
                % for categorical predictors betas need to be weighted by
                % normalized counts but since the predictors are expanded,
                % this count computation is the same as mean value in that
                % column, hence the same formula should hold, just that in
                % the end we need accumulate all the values belonging to a
                % single categorical predictor, that is done in the second
                % if-else block below

                % for all matrix data with categorical predictors,
                % variables need to be encoded prior to expanding
                if ~this.TableInput
                    thisX = classreg.learning.internal.encodeCategorical(thisX,this.BlackboxModel.VariableRange);
                    q = classreg.learning.internal.encodeCategorical(q,this.BlackboxModel.VariableRange);
                end
                thisX = classreg.learning.internal.expandCategorical(thisX,...
                    mdl.CategoricalPredictors,mdl.DataSummary.OrdinalPredictors,...
                    mdl.VariableRange);
                q = classreg.learning.internal.expandCategorical(q,...
                    mdl.CategoricalPredictors,mdl.DataSummary.OrdinalPredictors,...
                    mdl.VariableRange);
            end

            % Linear cshapley formula
            % For the interventional value function, linear cshapley values
            % for a predictor j reduce to the following expression
            %
            % phi(j) = beta*q(j) - beta*E(x(j))
            %
            % where q is the query point, and x is a random sample drawn from the
            % training distribution.
            % If x(j) is a continuous predictor, it is straightforward to
            % compute this expected value. It can be estimated to be the
            % average value taken by that predictor.
            % For a categorical predictor which is a discrete predictor the
            % first term is the evaluation of the query point which is the
            % scalar function value, and the second term is estimating the
            % expected value of a discrete variable which can take any of
            % the unique categorical values. This can be estimated as the
            % weighted average of the linear coefficients for each
            % category weighted by the normalized counts or frequency of
            % that category. Luckily, if we dummify the categorical
            % variable using one hot encoding the expression for the counts
            % matches the expression for the average or the mean of the
            % column. It is just that in the end the contribution of each
            % category needs to be accumulated to explain the contribution
            % of the categorical predictor as a whole. This step is
            % performed immediately after.
            if isprop(this.BlackboxModel, "Mu") && ~isempty(this.BlackboxModel.Mu)
                mu = this.BlackboxModel.Mu;
                sig = this.BlackboxModel.Sigma;
                q = semisupervised.Model.standardizePredictData(q,mu,sig);
                phi = coeff.*q; % mean of standardized X is zero
            else
                phi = coeff.*(q-mean(thisX)); % scores for the positive class
            end

            if hasCat
                % accumulate phi back for categorical predictors
                numcats = cellfun(@numel, mdl.VariableRange);
                indices = 1:M;
                numcats(~numcats) = 1; % account for the continuous predictor
                keys = repelem(indices,numcats);
                phi = accumarray(keys',phi')';
            end
            if this.Type == "classification"
                phi = [-phi; phi]; % same convention as classreg.learning.classif.CompactClassifByBinaryRegr.score
            end
            phi = phi';
            if returnObject
                this = makeShapleyValuesTable(this,phi);
            end
        end

        function [phi, this] = fitKernelSHAP(this,q,M,numsubsets, ...
                useparallel, doConditional,thisX, f, ...
                Xmat, qmat, tableinput, returnNaNShapleyValues, Z, wts)
            if nargin == 12 % when Z and wts are not precomputed, compute them
                [Z,wts] = cshapley.coalitionMatrix(M,numsubsets); % coalition matrix, Z and weights for each row in Z
            end
            returnObject = nargout==2;
            phi = nan(M, this.NumClasses); % initialize with nan
            if isempty(Z) % if no subsets besides the infinite weight subsets exist, return NaN
                if returnObject
                    this = makeNaNShapley(this);
                end
                return;
            end

            wts = wts./sum(wts); % normalize the weights
            if returnNaNShapleyValues
                % this implies one of the following
                % (1) the query point had nans as continuous predictors and distance computation was needed
                % (2) X has NaNs in all rows and the model cannot predict for NaN data i.e. non-tree models
                % return NaN Shapley values early
                if returnObject
                    this = makeNaNShapley(this);
                end
                return;
            end
            type = this.Type;
            if type =="classification"
                K = this.NumClasses;
            else % regression and function handle
                K = 1;
            end
            categcols = this.CategoricalPredictors;
            ev = cshapley.expectedValues(f, Z, type, K, thisX, q, Xmat, qmat, tableinput,categcols, useparallel, doConditional); % first column of Z is all ones (intercept)

            % account for the subsets of size M and 0 using variable
            % elimination and satisfying the sum constraint

            % 1. compute v0 and vM
            if type =="classification"
                [~,scoreQueryPoint] = f(q);
                v0 = this.Intercept; % the intercept
                vM = scoreQueryPoint;
            else
                v0 = this.Intercept; % the intercept
                vM = f(q);
            end

            % 2. Subtract the intercept from the ev
            v = ev - v0;
            fx = vM - v0;

            % 3. solve the least squares problem
            % min (Z*phi - v)*W*(Z*phi - v) s.t sum(phi) = fx
            v = v.*sqrt(wts); % make weighted v
            Z = Z.*sqrt(wts); % make weighted Z
            % explicitly remove the constraint
            zM = Z(:,end);
            Z = Z(:,1:end-1);
            Z = Z - zM;
            phiExceptLast = pinv(Z)*(v - zM.*fx);
            phi = [phiExceptLast;fx - sum(phiExceptLast,1)]; % do not forget to add 1 to the call in sum for the edge case where phiExceptLast is a row-vector (bivariate)
            if returnObject
                this = makeShapleyValuesTable(this,phi);
                % set the internal property UseParallel for debugging
                this.UseParallel = useparallel;
            end
        end

        function this = makeShapleyValuesTable(this,phi)
            % ======= make a table of cshapley values =========
            mdl = this.BlackboxModel;
            if ~this.IsFunctionHandle
                predictorNames = string(mdl.PredictorNames)';
            else
                if this.TableInput
                    predictorNames = string(this.X.Properties.VariableNames)'; % use table variable names as predictor names
                else
                    D = size(this.X,2);
                    predictorNames = string(classreg.learning.internal.defaultPredictorNames(D))';
                end
            end
            if this.Type=="classification" % ====== CLASSIFICATION ========
                classNames = strtrim(string(mdl.ClassNames))';
                predictorNamesTable = table(predictorNames, 'VariableNames', {'Predictor'}); % predictor
                phiTable = array2table(phi, 'VariableNames',classNames);
                this.ShapleyValues = [predictorNamesTable phiTable];
            else % ====== REGRESSION ========
                this.ShapleyValues = table(predictorNames,phi,'VariableNames',{'Predictor',...
                    'ShapleyValue'});
            end
        end

        function tbl = iChangeColumnNameForRegression(this, tbl)
            % change the column name from "ShapleyValue" to just "Value" to avoid
            % repetition
            if this.Type == "regression" && ~isempty(tbl)
                tbl = renamevars(tbl, "ShapleyValue", "Value");
            end
        end

        function [thisX, q, mdl, f, Xmat, qmat, tableinput, returnNaNShapleyValues, useKernelSHAP, skipShapleyCompute] = ...
                prepareData(this, mdl, q, doConditional, canUseTreeSHAP, hasSurrogateSplits,isCutVariable, isBatchFit)
            % This method is used for converting tabular data and removing
            % NaN data if necessary.
            %
            % thisX and q are inputs to scoring, these are the variables used for function handle evaluation or model evaluation
            % qmat and Xmat are vector/matrix forms of thisX and q
            %
            % These are used only when needed otherwise, they remain empty.
            % e.g. vanilla kernel shap with matrix data will not lead to creation of qmat and Xmat but for table inputs and vanilla kernel shap, Xmat and qmat are still needed for faster prediction
            %
            % However, conditional kernel shap needs some preprocessing on the data for distance computation.
            % Hence, for such cases, Xmat and qmat are always needed for modification before distance computation.
            if nargin == 4
                canUseTreeSHAP = false;
                isBatchFit = false;
                hasSurrogateSplits = false;
                isCutVariable = [];
            end
            N = size(q,1);
            skipShapleyCompute = false(N,1);
            returnNaNShapleyValues = false;
            useKernelSHAP = false;
            Xmat = [];
            qmat = [];
            isFunctionHandle = this.IsFunctionHandle;
            tableinput = this.TableInput;
            thisX = this.X(this.SampledObservationIndices, :);
            removedXNaNsAlready = false;
            categcols = this.CategoricalPredictors;
            % for tables, we must convert the original data to matrices for
            % the following purposes:
            % 1. Function handle input : if distance computation is needed i.e. the conditional approach
            % 2. Classreg model : in all cases, convert to matrices for speedy prediction

            isNonEmptyXmat = tableinput; % when Xmat is needed e.g. when the original input is a table

            if isFunctionHandle
                f = mdl;
            else
                f = @(x)predict(mdl,x);
            end

            istableq = istable(q);
            if tableinput % ============= TABLE X =========
                if isFunctionHandle
                    if ~doConditional % vanilla kernel shap
                        % Note that for vanilla kernel shap, Xmat need not be
                        % created unless q is not a table (then the function
                        % handle can evaluate on numeric matrices). For such
                        % cases convert thisX to a matrix for faster indexing.
                        % In all other cases for vanilla shap, the tables are
                        % directly passed to the evaluation.
                        if ~istable(q) % all-numeric q
                            Xmat = table2array(thisX); % table2array must be used because f can evaluate the numeric X as is and not with the encoded categories
                            % if q is a numeric matrix (and since f can predict on q, a numeric row)
                            thisX = Xmat; % this means X is an all-numeric table and can be used as a matrix to predict on matrix,
                            tableinput = false;
                        end
                        isNonEmptyXmat = false;
                    else % doConditional
                        % For the conditional method Xmat and qmat must be
                        % created for distance computation
                        %
                        % Note however, for function handles, we cannot use
                        % this  matrix for prediction, because if function
                        % handles were passed with tables, they are
                        % supposed to work with tables alone, not an
                        % encoding we use in the Statistics and Machine
                        % Learning Toolbox
                        [Xmat,~,vrange,~,args] = classreg.learning.internal.table2FitMatrix(thisX,[], 'CategoricalPredictors', categcols);
                        % For distance computation specifically we need
                        % to ensure that ordinals are not included
                        % (because no distance computation recipe
                        % exists), and throw an error for those
                        ords = args{6};
                        if any(ords)
                            error(message('stats:responsible:shapley:OrdinalsNotSupportedForConditional'));
                        end
                        prednames = args{4};

                        % the any(categcols) is not needed for function
                        % handles, it is just being added in the rarest of
                        % rare case where someone uses a classreg-like
                        % scoring function as a function handle and trains
                        % the model on a table, uses X as a table, and X
                        % has all numeric data but the user wishes to treat
                        % some of these numeric columns as categorical.
                        % Then this model will be able to predict on
                        % numeric matrices too (depends on the classreg
                        % function, but this is kind of an undocumented
                        % workflow, supported in classreg).
                        qmat = classreg.learning.internal.table2PredictMatrix(q,[],[],...
                            vrange,categcols, prednames,true);
                    end
                else %  A classreg model, we know the behavior of predict() and how to convert such tables
                    if doConditional
                        dataSummary = mdl.DataSummary;
                        if isfield(dataSummary,'OrdinalPredictors')
                            isOrdinal = any(dataSummary.OrdinalPredictors);
                        end
                        cshapley.checkForOrdinalsInData(thisX,isOrdinal); % this check is the same as the one done in fit() because the user can also modify "method" in fit()
                    end
                    % ===== MODEL TRAINED ON A MATRIX =====
                    if ~mdl.TableInput % remember this matrix could also have categoricals in it
                        % when trained on a matrix, and a valid table is
                        % passed i.e. a table with {'x1','x2',...} as its
                        % variable names,
                        % it must be an all-numeric table,
                        % convert it to a matrix for faster operations
                        Xmat = table2array(thisX(:,mdl.PredictorNames));
                        thisX = Xmat; % for faster prediction
                        if istableq
                            qmat = table2array(q(:,mdl.PredictorNames));
                            q = qmat; % for faster prediction
                        end
                        if doConditional % encode categories if computing distances
                            qmat = classreg.learning.internal.table2PredictMatrix(q,[],[],...
                                mdl.DataSummary.VariableRange,mdl.CategoricalPredictors, mdl.PredictorNames,true);
                            Xmat = classreg.learning.internal.table2PredictMatrix(thisX,[],[],...
                                mdl.DataSummary.VariableRange,mdl.CategoricalPredictors, mdl.PredictorNames,true);
                        end
                        tableinput = false;
                        if canUseTreeSHAP && ~doConditional % can directly use tree impl for trees instead of converting using fromStruct
                            useKernelSHAP = checkForSurrogates(hasSurrogateSplits, thisX, q,isCutVariable,isBatchFit);
                            if all(~useKernelSHAP) % if tree shap is the only algorithm needed return because can use the impl directly
                                return;
                            end
                        end
                    else % ===== MODEL TRAINED ON A TABLE =====
                        qmat = classreg.learning.internal.table2PredictMatrix(q,[],[],...
                            mdl.DataSummary.VariableRange,mdl.CategoricalPredictors, mdl.PredictorNames,true);
                        dataSummary = mdl.DataSummary;
                        Xmat = classreg.learning.internal.table2PredictMatrix(thisX,[],[],...
                            dataSummary.VariableRange,mdl.CategoricalPredictors, mdl.PredictorNames,true);
                        if canUseTreeSHAP && ~doConditional % can directly use tree impl for trees instead of converting using fromStruct
                            useKernelSHAP = checkForSurrogates(hasSurrogateSplits, Xmat, qmat,isCutVariable,isBatchFit);
                            if all(~useKernelSHAP) % if tree shap is the only algorithm needed return because can use the impl directly
                                thisX = Xmat; % use matrices for the builtin
                                q = qmat; % use matrices for the builtin
                                return;
                            end
                        end
                        if ~isempty(dataSummary.CategoricalPredictors) % convert model if there are categorical predictors
                            % for most tables, we should be able to convert both the table and the model itself
                            % since the model was trained on heterogeneous tables, the model needs to be modified to work with numeric matrices
                            % the following utility will help convert the model using toStruct/fromStruct
                            % this conversion is a MAJOR boost for performance for tabular data in some cases it is as high as 20X
                            [mdl, status] = cshapley.convertModelToPredictOnMatrix(mdl); % returns a compact model with a numeric variable range for categorical predictors
                            if status > 0 % model conversion was successful, use the converted matrices for everything
                                f = @(x)predict(mdl,x); % update the predict call to use the modified model
                                thisX = Xmat;
                                q = qmat;
                                tableinput = false;
                            else % for status == -1, the rare cases in which toStruct/fromStruct cannot help convert the model
                                % heterogeneous table, shrink it to just the predictor names
                                prednames = this.BlackboxModel.PredictorNames;
                                % shrink the tables to just predictor names
                                if istableq
                                    q = q(:,prednames);
                                    % if a numeric value is passed, it means data was an all-numeric table (since the model can predict on q)
                                end
                                thisX = thisX(:,prednames);
                            end
                        else
                            thisX = Xmat;
                            q = qmat;
                            tableinput = false;
                        end
                    end
                end
            else % ============ MATRIX X ===============
                if istableq % if X was not a table, and model can predict on q that means, q is an all numeric floating-type table with {'x1','x2', 'x3' .....}
                    if ~isFunctionHandle
                        qmat = table2array(q(:,mdl.PredictorNames));
                        q = qmat;
                    else
                        qmat = table2array(q); % when q is a table, we will need to create a copy into a vector qmat anyways
                        q = qmat;
                    end
                end
                if doConditional
                    isNonEmptyXmat = true;
                    % qmat is a vector which is used in distance
                    % computation
                    % Xmat is a matrix which is used in distance
                    % computation
                    % Ensure that their categories are encoded correctly
                    % for distance computation
                    if ~isFunctionHandle
                        % for classreg models with numeric matrices, use
                        % table2PredictMatrix to encode categories as
                        % integers, needed by internal.stats.heteropdist2
                        prednames = mdl.PredictorNames;
                        vrange = mdl.DataSummary.VariableRange;
                        Xmat = classreg.learning.internal.table2PredictMatrix(thisX,[],[],...
                            vrange,categcols, prednames,true);
                        qmat = classreg.learning.internal.table2PredictMatrix(q,[],[],...
                            vrange,categcols, prednames,true);
                    else
                        % We need not support the workflow for passing
                        % numeric matrices with categorical inputs for
                        % function handles
                        % We can always have the user pass tables for such
                        % use cases.
                        % However, just to be consistent with classreg, in
                        % case someone tries to pass purely numeric
                        % matrices with some columns to be treated as
                        % "categorical", do the same thing as classreg and
                        % encode the categories before distance
                        % computation for the "conditional" method.
                        % Of course, for such use cases, we are assuming
                        % that the function handle will be written in a
                        % way, that it can treat these columns as
                        % "categorical" in the function evaluation too i.e.
                        % the user is responsible for how the categories
                        % are handled in f, and we are responsible for how
                        % they are handled in the distance computation.
                        %
                        % To get the variable range, convert this matrix
                        % "thisX" to a table, and pass it to
                        % table2FitMatrix : this will populate the relevant
                        % arguments we can use immediately after to convert q
                        % Note that we are in the branch where thisX is a
                        % matrix, so we can call array2table on it.
                        [Xmat,~,vrange,~,args] = classreg.learning.internal.table2FitMatrix(array2table(thisX),[], 'CategoricalPredictors', categcols);
                        % For distance computation specifically we need
                        % to ensure that ordinals are not included
                        % (because no distance computation recipe
                        % exists), and throw an error for those
                        prednames = args{4};
                        qmat = classreg.learning.internal.table2PredictMatrix(q,[],[],...
                            vrange,categcols, prednames,true);
                    end
                end
                % check for surrogate splits in the case of a matrix as
                % well
                if canUseTreeSHAP && ~doConditional % can directly use tree impl for trees instead of converting using fromStruct
                    useKernelSHAP = checkForSurrogates(hasSurrogateSplits, thisX, q,isCutVariable,isBatchFit);
                    if all(~useKernelSHAP) % if tree shap is the only algorithm needed return because can use the impl directly
                        return;
                    end
                end
            end
            % ===== score related NaN checks ====
            % for non-trees, remove rows with any NaNs in them
            % this is because those models do not know how to
            % predict with missing data and would produce NaN
            % values in the expected value computation
            % (remove NaN values for function handles too)
            if isFunctionHandle || (~isFunctionHandle && ~cshapley.isTreeEnsembleGAM(mdl))
                if isFunctionHandle % for function handles, remove NaN rows ONLY if they produce NaN scores
                    scoresWithPossibleNaNs = f(thisX);
                    nanrows = isnan(scoresWithPossibleNaNs);
                    if isNonEmptyXmat % for function handles, there might still exist continuous predictors that are NaN and distance computation is needed (which should be later removed
                        Xmat(nanrows,:) = [];
                    end
                    thisX(nanrows,:) = [];
                    removedXNaNsAlready = false;
                else % for non-tree models in the toolbox e.g. SVM
                    if isNonEmptyXmat
                        nans = isnan(Xmat);
                    else
                        nans = isnan(thisX);
                    end
                    nanrows = any(nans,2);
                    if isNonEmptyXmat
                        Xmat(nanrows,:) = [];
                    end
                    thisX(nanrows,:) = [];
                    removedXNaNsAlready = true;
                end
            end
            if isempty(thisX)
                returnNaNShapleyValues = true;% all rows have nans, return NaN cshapley values
                skipShapleyCompute = true(N,1);
            end

            % If the query point has missing values for specific regression
            % models, return all NaNs
            if ~isFunctionHandle
                qmat = classreg.learning.internal.table2PredictMatrix(...
                    q,[],[],mdl.DataSummary.VariableRange,...
                    mdl.DataSummary.CategoricalPredictors,mdl.PredictorNames,true);
                qnan = isnan(qmat);
                hasNaNColumns = any(qnan,2);
                if any(hasNaNColumns)
                    modelsWithPredictionForMissing = classreg.learning.internal.modelsAllowPredictionForMissing;
                    for j=1:numel(modelsWithPredictionForMissing)
                        if isa(mdl,modelsWithPredictionForMissing{j})
                            returnNaNShapleyValues = true;
                            skipShapleyCompute = hasNaNColumns;
                        end
                    end
                end
            end

            % ===== distance calculation related NaN checks ====
            %
            % Xmat and qmat are used in distance calculation for the
            % conditional approach,
            %
            % only internal.stats.heteropdist2 knows how to handle nan
            % categories for categorical predictors,
            %
            % for any continuous predictor, nan rows must be removed
            % because continuous predictor distance calculation functions
            % (pdist2 and heteropdist2, continuous part) do not know how to
            % handle NaNs.
            %
            if doConditional
                if isempty(qmat) % if a numeric vector was passed for q
                    qnan = isnan(q);
                else
                    qnan = isnan(qmat); % for function handles we do not know if a nan query point actually results in a nan score, so we cannot say that we already checked for it
                end
                M = numPredictors(this);
                contidx = setdiff(1:M,categcols);
                if ~isempty(contidx) % remove rows with NaN continuous predictors
                    hasNaNContinuousColumns = any(qnan(:,contidx),2);
                    if doConditional && any(hasNaNContinuousColumns)
                        % in this case the score is non-NaN for a query point with NaN continuous predictors, for such query points we cannot do distance calculation
                        returnNaNShapleyValues = true;% all rows have nans, return NaN cshapley values
                        skipShapleyCompute = hasNaNContinuousColumns;
                    end
                    if ~removedXNaNsAlready % for tree based models, still remove continuous predictors that have NaNs
                        if isNonEmptyXmat
                            nans = isnan(Xmat);
                        else
                            nans = isnan(thisX); % X is a matrix
                        end
                        nans = nans(:,contidx);
                        nanrows = sum(nans,2)~=0;
                        if isNonEmptyXmat
                            Xmat(nanrows,:) = [];
                        end
                        thisX(nanrows,:) = [];
                        if isempty(thisX)
                            returnNaNShapleyValues = true; % all rows have nans, return NaN cshapley values
                            skipShapleyCompute = true(N,1);
                        end
                    end
                end
                % for all continuous data, standardize the data
                % for distance computation: (standardized euclidean distance)
                if all(~categcols)
                    M = numPredictors(this);
                    isCategorical = false(1,M);
                    if isNonEmptyXmat
                        [Xmat,mu,sigma] = semisupervised.Model.standardizeData(Xmat,isCategorical); % if the original data is a table
                    else
                        [Xmat,mu,sigma] = semisupervised.Model.standardizeData(thisX,isCategorical); % if the original data is a matrix
                    end
                    if ~isempty(qmat)
                        qmat = semisupervised.Model.standardizePredictData(qmat,mu,sigma);
                    else
                        qmat = semisupervised.Model.standardizePredictData(q,mu,sigma);
                    end
                end
            end
        end

        function [method,numsubsets,useparallel,M,istree,isens,islinear,hasSurrogateSplits,isCutVariable,userProvidedUseParallel] = parseAndValidateFitArgs(this, varargin)
            % ========== parse all the fit args ==========
            thisFileName = mfilename;
            fitparser = inputParser;
            % validate types
            validateStringScalar= @(x) validateattributes(x,{'string', 'char'}, ...
                {'scalartext','nonempty'});
            if(this.Method == "conditional-kernel")
                defaultMethod = 'conditional';
            else
                defaultMethod = 'interventional';
            end
            fitparser.addParameter('Method', defaultMethod ,validateStringScalar);
            fitparser.addParameter('MaxNumSubsets', this.NumSubsets, @(x)validateattributes(x,{'numeric'},...
                {'scalar','>',1,'finite', 'integer'}));
            fitparser.addParameter('UseParallel', false, @(x)validateattributes(x,{'logical','double','string','char'},{'vector'}));
            M = numPredictors(this);
            % ========== parse all the fit args ==========
            fitparser.parse(varargin{:});
            parsedResults = fitparser.Results;
            method = parsedResults.Method;
            numsubsets = min(parsedResults.MaxNumSubsets,2^M);
            useparallel = parsedResults.UseParallel;
            % ======validate some of these inputs ========
            % 1.method
            oldValidMethods = ["interventional-kernel", "conditional-kernel"]; % before 23a
            validMethods = ["interventional", "conditional"];
            if ~ismember('Method', fitparser.UsingDefaults) % user passed method
                if ismember(method, oldValidMethods)
                    error(message('stats:responsible:shapley:UseInterventionalOrConditional'));
                end
                method = validatestring(method,validMethods, thisFileName);
            end
            % tree specific outputs
            hasSurrogateSplits = false;
            isCutVariable = [];
               % 2. allow 0 and 1 as UseParallel
            useparallel = internal.stats.parseOnOff(useparallel,'UseParallel');
            userProvidedUseParallel = ~ismember('UseParallel', fitparser.UsingDefaults);
            if method == "interventional"
                % add warnings to say that "MaxNumSubsets"/"UseParallel"
                % is a no-op for some models
                learner = this.BlackboxModel;
                islinear = cshapley.isCandidateForLinearSHAP(learner);
                [istree, isens, hasSurrogateSplits, isCutVariable] = cshapley.isCandidateForTreeSHAP(learner, numPredictors(this));

                if istree || isens || islinear
                    if ~ismember('MaxNumSubsets', fitparser.UsingDefaults) % MaxNumSubsets is a hyperparameter for kernelSHAP, use kernelSHAP even when specified for tree ensembles, learners and linear learners
                        if islinear
                            warning(message('stats:responsible:shapley:MaxNumSubsetsProvidedForLinear'));
                            islinear = false;
                        else % istree and isens
                            warning(message('stats:responsible:shapley:MaxNumSubsetsProvidedForTrees'));
                            istree = false;
                            isens = false;
                        end
                    end
                     
                    if this.IsLocal % 24a: for multiple query points, we do support this parameter
                        if userProvidedUseParallel
                            warning(message('stats:responsible:shapley:UseParallelNotApplicableToTreesLinear'));
                        end
                    end
                else
                    % the expected value computation for kernel methods (the most expensive step) has
                    % a complexity of computing size(thisX,1)*NumSubsets
                    % predictions, this scales up pretty fast e.g. for a 4
                    % predictor X with a 1000 observations, this amounts to 16000
                    % predictions.
                    % suggest the user either downsampling or using the
                    % 'UseParallel' flag
                    if this.NumObservationsToSample > 1000
                        if ~useparallel
                            warning(message('stats:responsible:shapley:BigXSlowComputation'));
                        else
                            warning(message('stats:responsible:shapley:BigXSlowComputationSuggestSubsample'));
                        end
                    end
                end
            else
                istree = false;
                isens = false;
                islinear = false;
            end
        end 

        function M = numPredictors(this)
            if ~this.IsFunctionHandle
                M = numel(this.BlackboxModel.PredictorNames); % M is the number of features, N is the number of observations
            else
                M = size(this.X,2);
            end
        end

        function this = makeNaNShapley(this)
            M = numPredictors(this);
            if this.Type == "classification"
                phi = nan(M,this.NumClasses);
                this = makeShapleyValuesTable(this,phi);
            else
                phi = nan(M,1);
                this = makeShapleyValuesTable(this,phi);
            end
        end 


        function barObjects = localPlot(this,ax,isRegression,queryPointIndices,M,K,m,cnames,classidx,userPassedClassNames,userClassNames)
            if this.IsLocal
                offset = 1; % first column is features
                phi = this.ShapleyValues{:,classidx+offset};
            else
                phi = this.RawShapleyValues(queryPointIndices,:,classidx);
                phi = reshape(phi,M,K);
            end
            absolutePhi = abs(phi);

            if isscalar(classidx) % when there is only one column of phis to be plotted
                [~,predictorOrder] = sort(absolutePhi,'descend');
            else
                sums = sum(absolutePhi,2);
                [~,predictorOrder] = sort(sums,'descend');
            end
            predidx = predictorOrder(1:m);
            % bar plots from the bottom, flip the values, so that the top
            % values appear at the top
            phiToPlot = flipud(phi(predidx,:));
            preds = flipud(this.ShapleyValues.Predictor(predidx)); % top m predictors
            barObjects = barh(ax,phiToPlot);
            yticks(ax,1:m);
            yticklabels(ax,preds);

            % add title, labels, subtitle and if necessary, legend
            title(ax,getString(message('stats:responsible:shapley:PlotTitle')));
            set(ax,'TickLabelInterpreter', 'none'); % set the tick label interpreter to none
            xlabel(ax,getString(message('stats:responsible:shapley:XLabel')));
            ylabel(ax,getString(message('stats:responsible:shapley:YLabel')));
            if isRegression
                if this.IsLocal
                    queryPointPrediction = num2str(this.BlackboxFitted); % for regression, this is the response
                else
                    queryPointPrediction =  num2str(this.BlackboxFitted(queryPointIndices));
                end
                averagePrediction = num2str(this.Intercept);
                subtitle(ax, sprintf(getString(message('stats:responsible:shapley:PlotSubtitleForRegression',  queryPointPrediction, averagePrediction))), 'Interpreter','none');
            else % for classification use the score of the predicted class
                if this.IsLocal
                    predictedClass = this.BlackboxFitted;
                else
                    predictedClass = this.BlackboxFitted(queryPointIndices);
                end
                [~,predictedClassIdx] = ismember(classreg.learning.internal.ClassLabel(predictedClass), cnames);
                [~,scores] = predict(this.BlackboxModel, this.QueryPoints);
                queryPointPrediction = num2str(scores(queryPointIndices,predictedClassIdx));
                averagePrediction = num2str(this.Intercept(:, predictedClassIdx));
                predictedClassString = string(predictedClass);
                subtitle(ax,sprintf(getString(message('stats:responsible:shapley:PlotSubtitleForClassification', predictedClassString,  queryPointPrediction, averagePrediction))),'Interpreter','none');
                if userPassedClassNames % include the legend only when classnames are passed
                    l = legend(ax,string(labels(userClassNames)),'Location','best','Interpreter','none'); % set the legend interpreter to none
                    title(l,getString(message('stats:responsible:shapley:LegendTitle')));
                end
            end
        end
    end

    methods % get-set methods
        function arr = get.RawShapleyValues(this)
            rawtbl = this.ShapleyValues;
            N = size(this.QueryPoints,1);
            M = numPredictors(this);
            K = this.NumClasses;
            arr = nan(N,M,K,'like',rawtbl{1,2}); % use one of the values for the "type"
            offset = 1; % first column in the table is predictors
            for k = 1:K
                arr(:,:,k) = rawtbl{:,k+offset}';
            end
        end

        function tbl = get.Shapley(this)
            tbl = this.ShapleyValues;
            tbl = iChangeColumnNameForRegression(this, tbl);
        end
        
        function this = set.QueryPoint(this,q)
            this.QueryPoints = q;
        end

        function q = get.QueryPoint(this)
            q = this.QueryPoints;
        end

        function islocal = get.IsLocal(this)
            islocal = size(this.QueryPoints,1) == 1;
        end

        function tbl = get.MeanAbsoluteShapley(this)
            rawtbl = this.ShapleyValues;
            if ~isempty(rawtbl)
                varnames = rawtbl.Properties.VariableNames;
                tbl = table();
                tbl.(varnames{1}) = rawtbl.(varnames{1}); % first column
                % cannot use varfun because it changes the predictor order
                K = this.NumClasses;
                offset = 1;
                for k = 1:K
                    classidx = offset+k;
                    shapvals = rawtbl{:,classidx};
                    meanabsshap = mean(abs(shapvals),2,"omitnan");
                    tbl.(varnames{classidx}) = meanabsshap;
                end
            else
                tbl = [];
            end
            tbl = iChangeColumnNameForRegression(this, tbl);
        end
    end

    % helper methods: use static, hidden methods to keep these utilities
    % contained in this file (easy to maintain) and be tested independently
    methods (Static, Hidden)
        function ev = expectedValues(f, Z, type,K, thisX, q, Xmat, qmat, tableinput,categcols, useparallel, doConditional)
            % EXPECTEDVALUES returns the expected value of the scores for
            % all subsets The MOST EXPENSIVE STEP BY FAR in computing
            % Shapley Values (99% of the time is spent here for most
            % problems).
            %
            % ev is the expected value of the prediction of artificial
            % samples created for the query point. These artificial samples
            % are created in such a way that if a feature is included in a
            % subset => it is taken from the query point if a feature is
            % not included in a subset, it is taken from the training set
            %
            % f is the scoring/response function
            % Z is the coalition matrix
            % type is either 'classification' or 'regression'
            % K is the number of classes/columns in ShapleyValues
            % thisX is the value of X stored on the object
            % q is the query point
            % Xmat is the output matrix from table conversion of X
            % qmat is the output vector from table conversion of q
            % tableinput is true if the user originally passed X as a table
            % categcols is a vector of indices of categorical predictors
            % useparallel is a logical scalar that determines if parfor be used over for
            % doConditional is a logical scalar that determines if the conditional expected values are computed
            if isempty(Xmat)
                isNonEmptyXmat = false;
            else
                isNonEmptyXmat = true;
            end
            numsubsets = size(Z,1);
            N = size(thisX,1);

            % parfor needs the following variable to be declared, even for the non-conditional
            % approach
            numNeighbors = ceil(0.1*N);  % consider 10% of the neighbors for the conditional approach
            allContinuous  = all(~categcols);

            % ==== A NOTE on the Parallelization Scheme ====
            %
            % Though the total number of atomic operations involved are the
            % same, i.e we need to predict on numsubsets*numobservations
            % samples, it is still known that predict for classreg
            % functions can predict on multiple samples at once, we must
            % try to minimize the number of predict() calls made in our
            % loops.
            %
            % Profiling reveals that predict calls in the for loop are the
            % most expensive operation, and reducing the number of these
            % predict calls shows better performance. For example, a
            % problem with 4 features will have a coalition matrix of size
            % (2^4) i.e. numsubsets would be 16, and yet it can have any
            % number of observations, likely to be much higher than 16, say
            % 2000, for such cases it makes sense to predict on a matrix of
            % size 2000-by-4 inside the loop.
            %
            % On the other hand, we had 50 features, the exact problem of
            % size 2^50 is not feasible to solve, but the numsubsets
            % allocated should still be high to get good accuracy on
            % cshapley values and let's say the numsubsets is about 10,000,
            % and the number of observations is much fewer say 2000. For
            % this case, we should predict on a matrix of size 10000-by-50.
            %
            % Note that this strategy benefits both serial and parallel
            % versions: doing more work for each iteration makes each
            % worker do more work in parallel, and for serial, it reduces
            % the calling overhead due to multiple predict() calls.

            loopOverSubsets = tableinput || (numsubsets < N) || doConditional ; % a boolean that is set to true when the strategy of looping over subsets is used
            % it is false when looping over observations, only used for vanilla kernel shap using small matrices of observations with many subsets
            if ~useparallel % ======== SERIAL IMPLEMENTATION ========
                if loopOverSubsets %  ========== table input, small budget many observations ============
                    meanScore = zeros(numsubsets,K);
                    for coalition = 1:numsubsets
                        meanScore(coalition,:) = cshapley.meanScoreForLoopOverSubsets(f,coalition,N, type,Z,thisX,q,isNonEmptyXmat, Xmat,qmat,categcols,doConditional, numNeighbors, allContinuous);
                    end
                    ev = meanScore;
                else  % ========== matrix input, big budget small number of observations ============
                    Q = repmat(q,numsubsets,1);
                    QZ = Q(Z);
                    sumScoresForAllCombinations = zeros(numsubsets,K);
                    for obsidx = 1:N
                        sumScoresForAllCombinations = sumScoresForAllCombinations + cshapley.scoresForLoopOverObservations(f,obsidx,numsubsets,type, thisX, QZ,Z);
                    end
                    ev = sumScoresForAllCombinations/N; % take the average
                end

            else % ======== PARALLEL IMPLEMENTATION ==========
                if loopOverSubsets %  ========== table input, small budget many observations ============
                    meanScore = zeros(numsubsets,K);
                    parfor coalition = 1:numsubsets
                        meanScore(coalition,:) = cshapley.meanScoreForLoopOverSubsets(f,coalition,N, type,Z,thisX,q,isNonEmptyXmat, Xmat,qmat,categcols,doConditional, numNeighbors, allContinuous);
                    end
                    ev = meanScore;
                else  % ========== matrix input, big budget small number of observations ============
                    Q = repmat(q,numsubsets,1);
                    QZ = Q(Z);
                    sumScoresForAllCombinations = zeros(numsubsets,K);
                    parfor obsidx = 1:N
                        sumScoresForAllCombinations = sumScoresForAllCombinations + cshapley.scoresForLoopOverObservations(f,obsidx,numsubsets,type, thisX, QZ,Z);
                    end
                    ev = sumScoresForAllCombinations/N; % take the average
                end
            end
        end

        function meanScore = meanScoreForLoopOverSubsets(f,coalition,N, type,Z,thisX,q,isNonEmptyXmat, Xmat,qmat,categcols,doConditional,numNeighbors, allContinuous)
            % LOOP BODY FOR LOOPING OVER SUBSETS
            % Useful for the following cases:
            % 1. matrix with large number of observations and small number
            %    of subsets
            % 2. the conditional approach
            % 3. table inputs
            z = Z(coalition,:);
            % where a predictor is chosen, take its value from Q => Q(z) must survive
            % where it is not chosen, take its value from X => X(~z) must survive
            Xprime = thisX;
            Xprime(:,z) = repmat(q(:,z), N, 1);
            if doConditional
                if allContinuous % purely continuous data, already standardized, in prepareData, Xmat will always exist for all continuous predictors because it will be the standardized values in X
                    distances = pdist2(Xmat(:,z),qmat(z));
                else % heterogeneous data
                    % scan the subspace to see if there are any
                    % categorical columns in the subspace, and
                    % pass those indices to heterpdist2
                    subspaceidx = find(z);
                    [~,categoricalsForSubspace] = ismember(categcols,subspaceidx);
                    categoricalsForSubspace = categoricalsForSubspace(categoricalsForSubspace~=0);
                    if isNonEmptyXmat
                        distances = internal.stats.heteropdist2(Xmat(:,z),qmat(z), 'goodall3',categoricalsForSubspace); % use goodall for speed
                    else % will be used for matrices with numeric categorical predictors
                        distances = internal.stats.heteropdist2(thisX(:,z),q(z), 'goodall3',categoricalsForSubspace); % use goodall for speed
                    end
                end
                [~, neighboridx] = mink(distances, numNeighbors);
                Xprime = Xprime(neighboridx,:); % include only neighbors for scoring for the conditional approach
            end
            % not conditional, vanilla kernel shap
            if type=="classification" %  ========== CLASSIFICATION ============
                [~, scores] = f(Xprime);
            else % ========== REGRESSION ============
                scores = f(Xprime);
            end
            meanScore = mean(scores,1); % add this one for the edge case of X being of size 1 (mean will shrink it to a scalar)
        end

        function scores = scoresForLoopOverObservations(f,obsidx,numsubsets,type, thisX, QZ,Z)
            % LOOP BODY FOR LOOPING OVER OBSERVATIONS
            % Useful only for the following case:
            % matrix input with small number of observations and a large
            % number of numsubsets (wide data)
            Xprime = repmat(thisX(obsidx,:),numsubsets,1);
            % where a predictor is chosen, take its value from Q => Q(Z) must survive
            % where it is not chosen, take its value from X => X(~Z) must survive
            Xprime(Z) = QZ;
            if type=="classification" %  ========== CLASSIFICATION ==========
                [~, scores] = f(Xprime);
            else % ========== REGRESSION ============
                scores = f(Xprime);
            end
        end


        function [Z, subsetWeights] = coalitionMatrix(M, numsubsets)
            % The purpose of this function is to return an appropriately sized
            % coalition matrix of trues and falses (0s and 1s).
            %
            % The performance of this function is NOT CRITICAL at all (99% of
            % the time is spent in the expected value computation). Hence, the
            % UseParallel flag has no implications for this helper.
            %
            % Z is a NumSubsets-by-M coalition matrix, where M is the
            % actual number of DIMENSIONS in the problem numsubsets is the number of
            % subsets to consider for this problem.
            %
            % The NumSubsets property can also be interpreted as the level of
            % approximation. The solution for an exact Shapley value
            % computation will require Z to be of size 2^M-by-M.

            if numsubsets ==2  % numsubsets is always greater than 1 (already validated in fit())
                % there are two subsets with infinite weights, they will be added to Z before solving the least squares
                % return early for this case of a budget of only two coalitions
                Z = [];
                subsetWeights = [];
                warning(message('stats:responsible:shapley:MaxNumSubsetsTooSmall'));
                return;
            end

            if M==2 % only two features, spell out Z
                switch numsubsets
                    case 3
                        Z = [true false];
                        subsetWeights = 1/2;
                    case 4
                        Z = [true false
                            false true];
                        subsetWeights = [1/2; 1/2];
                end
            else
                numsubsetsLargeEnough = false; % numsubsets is large enough to enumerate all subsets (2^M - 2 subsets)
                hasOddNumDimensions = mod(M,2)==1; % for odd number of features every single coalition has a complement e.g. (3choose0) + (3choose1) | (3choose2) + (3chose3)
                halfNumDimensions = floor(M/2); % symmetric problem
                featureSubsetSizes = uint64(1:halfNumDimensions)'; % a column of size halfNumDimensions
                binomialCoefficientsHalf= arrayfun(@(k)nchoosek(M,k), featureSubsetSizes); % a column of size halfNumDimensions
                INFINITY = intmax("uint64");
                cumulativeEnumerationsHalf = cumsum(binomialCoefficientsHalf); % same size as halfNumDimensions,
                % can contain overflown integers but that is ok, as long as the budget is reasonable, but if the specified budget also exceeds intmax, throw an error
                if any(cumulativeEnumerationsHalf) == INFINITY && budget > INFINITY
                    error(message('stats:responsible:shapley:InfiniteNumSubsets'));
                end
                % ideally we want to enumerate all the subsets,
                totalSubsetsWhole = 2^M; % this is the total number of subsets
                if numsubsets >= totalSubsetsWhole  % if the numsubsets exceeds the total number of enumerations to solve this problem, we can obtain an exact solution
                    numsubsetsLargeEnough = true;
                else % numsubsets is less than required for enumeration of all possible subsets
                    if numsubsets == 2*M + 2 % the case when subsets of size 1 can be enumerated
                        featureSubsetSizeExact = 1;
                    elseif numsubsets > 2*M + 2 % the case when more subsets than just size 1 can be enumerated
                        for ind = 1:halfNumDimensions
                            if numsubsets < 2*cumulativeEnumerationsHalf(ind) + 2
                                featureSubsetSizeExact = featureSubsetSizes(ind-1); % number of pairs of subsets we can fully enumerate within the numsubsets
                                break;
                            end
                        end
                    else % for numsubsets less than 2*M + 2, here not all subsets of size 1 can be enumerated
                        numSubsetsMinus2 = numsubsets-2;
                        Z = false(numSubsetsMinus2,M);
                        subsetWeights = zeros(numSubsetsMinus2,1);
                        if numSubsetsMinus2 <= M % numsubsets are in the range (2, M+2]
                            Z(1:numSubsetsMinus2+1:numSubsetsMinus2*numSubsetsMinus2) = true;
                            subsetWeights(1:numSubsetsMinus2) = 1/M;
                        else % numsubsets in the range (M+2, 2M+2)
                            Z(1:numSubsetsMinus2+1:end) = true;
                            remainingSubsets = numSubsetsMinus2 - M;
                            Z(M+1:M+remainingSubsets,:) = ~Z(1:remainingSubsets,:); % fill the remaining with as many complements as you can
                            subsetWeights = (1/M)*(ones(numSubsetsMinus2,1));
                        end
                        warning(message('stats:responsible:shapley:MaxNumSubsetsTooSmall'));
                        return;
                    end
                end

                kernelWeightsHalf = zeros(halfNumDimensions,1);
                binomialCoefficientsHalf = double(binomialCoefficientsHalf); % cast to double for weight calculation
                for ind = 1:halfNumDimensions % here ind is |z| and allPossibleSubsetSizesHalf is values of (M choose |z|) computed apriori, essentially make a hash table of kernel weights
                    kernelWeightsHalf(ind)= (M-1)/(binomialCoefficientsHalf(ind)*ind*(M-ind)); % (M-1)/((M choose |z|)*(|z|)*(M-|z|), where is z is the subset, and |z| is the subset size
                end

                if numsubsetsLargeEnough % we can fully enumerate all the subsets if the numsubsets is large enough to accommodate the size
                    exactSize = 2^M-2; % there two subsets which have all features included and none of the features included, handle them separately
                    % coalition matrix, Z (exactSize-by-M)
                    %   
                    % +----------------------------------------+
                    % |        subsets of size 1               |
                    % |                                        |
                    % +----------------------------------------+
                    % |        subsets of size 2               |
                    % |                                        |
                    % +----------------------------------------+
                    % |                 .                      |
                    % |                 .                      |
                    % |       goes on exactSize/2 times        |
                    % |                 .                      |
                    % |                 .                      |
                    % +----------------------------------------+
                    % |        subsets of size M-1             |
                    % |       (complements of above)           |
                    % +----------------------------------------+
                    % |         subsets of size M-2            |
                    % |         (complements of above)         |
                    % +----------------------------------------+
                    % |                 .                      |
                    % |                 .                      |
                    % |       goes on exactSize/2 times        |
                    % |                 .                      |
                    % |                 .                      |
                    % +----------------------------------------+
                    Z = false(exactSize, M);
                    subsetWeights = zeros(exactSize,1);
                    % work on subsets of size 1 through ceil(M/2) (illustrated in the top half of the picture above)
                    % the rest will be complements of what you populate

                    % as written, every iteration of the following for-loop is independent i.e. we can obtain these enumerations for a given subset size independently
                    % however, parfor has very restrictive indexing, so the following for-loop cannot be changed to "parfor" right away
                    % a parallel approach to this problem needs to be thought out

                    for featureSubsetSize = 1:halfNumDimensions
                        hotIndices = nchoosek(1:M, featureSubsetSize); % hot indices are of size n-by-M where n = M-choose-featureSubsetSize
                        rowSize = binomialCoefficientsHalf(featureSubsetSize);
                        if featureSubsetSize ==1
                            rowOffset = 0; % avoid zero indexing into cumulativeEnumerations
                        else
                            rowOffset = cumulativeEnumerationsHalf(featureSubsetSize-1); % totalPossibleEnumerations are the cumulative sum so far
                        end

                        kw = kernelWeightsHalf(featureSubsetSize);
                        for row = 1:rowSize
                            Z(row+rowOffset, hotIndices(row,:)) = true; % flip the false bits in this subset
                            subsetWeights(row+rowOffset) = kw; % kernel weights are symmetric so fill all of them out
                        end
                    end
                    % for the subsets of size M-1 through M-ceil(M/2), do not do any work except taking the complement (bottom half of the above picture)
                    if hasOddNumDimensions % for odd number of dimensions every subset has its complement
                        Z(exactSize/2+1:end,:) = ~Z(1:exactSize/2,:); % add the complements to Z
                        subsetWeights(exactSize/2+1:end) = subsetWeights(1:exactSize/2); % the weights are the same
                    else
                        % for even number of dimensions, subsets of size ceil(M/2) do not have a complement
                        % the interval of these asymmetric indices is as follows:
                        rowEndSymmetryIdx = cumulativeEnumerationsHalf(end-1);
                        rowBeginSymmetryAgainIdx = cumulativeEnumerationsHalf(end)+1;
                        Z(rowBeginSymmetryAgainIdx:end,:) = ~Z(1:rowEndSymmetryIdx,:); % add the complements to Z
                        subsetWeights(rowBeginSymmetryAgainIdx:end) = subsetWeights(1:rowEndSymmetryIdx); % add the complements to Z
                    end

                else % do budgeted shap, remember that we can enumerate featureSubsetSizeExact fully but the rest can only be partially enumerated

                    % Use a deterministic approach to make the coalition matrix, we know
                    % that the highest weighted subsets are the most important, so try to
                    % put as many highest weighted subsets in the coalition matrix as we
                    % can. There is no random sampling going on here for this reason. The
                    % highest weighted subsets contribute most to the cshapley sum. The SHAP
                    % package uses a combination of this deterministic approach and random
                    % sampling for the last remaining subsets. It is just faster to not
                    % randomly sample even for the last few remaining subsets.

                    % coalition matrix, Z (numsubsets-by-M),
                    %
                    % +----------------------------------------+
                    % |        subsets of size 1               |
                    % |                                        |
                    % +----------------------------------------+
                    % |        subsets of size 2               |
                    % |                                        |
                    % +----------------------------------------+
                    % |                 .                      |
                    % |                 .                      |
                    % | goes on featureSubsetSizeExact/2 times |
                    % |                 .                      |
                    % |                 .                      |
                    % +----------------------------------------+
                    % |--   full subset enumerations above   --|
                    % +----------------------------------------+
                    % |                                        |
                    % |leftover numsubsets: partial enumeration|
                    % |                                        |
                    % |                                        |
                    % +----------------------------------------+
                    % |--       complements below            --|
                    % +----------------------------------------+
                    % |        subsets of size M-1             |
                    % |       (complements of above)           |
                    % +----------------------------------------+
                    % |         subsets of size M-2            |
                    % |         (complements of above)         |
                    % +----------------------------------------+
                    % |                 .                      |
                    % |                 .                      |
                    % | goes on featureSubsetSizeExact/2 times |
                    % |                 .                      |
                    % |                 .                      |
                    % +----------------------------------------+
                    numSubsetsMinus2 = numsubsets-2; % exclude the two trivial subsets (taking all features or taking none), will be added later by fitKernelSHAP
                    Z = false(numSubsetsMinus2, M);
                    subsetWeights = zeros(numSubsetsMinus2,1);

                    for featureSubsetSize = 1:featureSubsetSizeExact
                        hotIndices = nchoosek(1:M, featureSubsetSize); % hot indices are of size n-by-M where n = M-choose-featureSubsetSize
                        rowSize = binomialCoefficientsHalf(featureSubsetSize);
                        if featureSubsetSize ==1
                            rowOffset = 0; % avoid zero indexing into cumulativeEnu
                        else
                            rowOffset = cumulativeEnumerationsHalf(featureSubsetSize-1); % totalPossibleEnumerations are the cumulative sum so far
                        end
                        kw = kernelWeightsHalf(featureSubsetSize);
                        for row = 1:rowSize
                            Z(row+rowOffset, hotIndices(row,:)) = true; % flip the false bits in this subset
                            subsetWeights(row+rowOffset) = kw; % kernel weights are symmetric so fill all of them out
                        end
                    end

                    % now do partial enumeration for the next subset in line that would have been fully enumerated if we had more numsubsets
                    halfBudget = ceil(numSubsetsMinus2/2);
                    rowEndSymmetryIdx = floor(numSubsetsMinus2/2); % take complements of subsets up to here
                    rowBeginSymmetryAgainIdx = halfBudget+1; % populate complements starting at this index
                    partialEnumerationBlock = row+rowOffset+1:halfBudget;
                    hotIndices = cshapley.partiallyEnumeratedSubsets(1:M,featureSubsetSizeExact+1,numel(partialEnumerationBlock)); % this is the subset next in line to enumerate

                    for row = 1:numel(partialEnumerationBlock)
                        Z(partialEnumerationBlock(row), hotIndices(row,:)) = true;
                        subsetWeights(partialEnumerationBlock(row)) = kernelWeightsHalf(featureSubsetSizeExact+1);
                    end

                    Z(rowBeginSymmetryAgainIdx:end,:) = ~Z(1:rowEndSymmetryIdx,:); % add the complements to Z
                    subsetWeights(rowBeginSymmetryAgainIdx:end) = subsetWeights(1:rowEndSymmetryIdx); % add the complements to Z
                end
            end
        end

        function P = partiallyEnumeratedSubsets(predictoridx,subsetSize, numCombinations)
            %   partiallyEnumeratedSubset returns only up to numCombinations combinations instead of all combinations
            % This helper is used only to fill out the partialEnumerationBlock in the coalition matrix
            % For example, let's say we have 10 features and MaxNumSubsets is chosen to be 25
            % We know that with that budget, we can fully enumerate 22 subsets
            % i.e. the highest weighted subsets (of size zero, and complements, of size 1 and complements):
            % include all predictors, exclude all predictors, include only one predictor, exclude only one predictor
            % For the remaining 3 subsets, we will enumerate the NEXT subset size (subsets of size 2 and complements)
            % For this example the NEXT subset size, (featureSubsetSizeExact+1), is 2.
            % Since the partially enumerated size will be a subset of this known size, kernelWeightsHalf(featureSubsetSizeExact+1)

            predictoridx = predictoridx(:).'; % Make sure it is a row vector
            n = length(predictoridx);

            P = zeros(numCombinations, subsetSize, 'like', predictoridx);

            % Compute P one row at a time:
            ind = 1:subsetSize;
            P(1, :) = predictoridx(1:subsetSize);
            for i=2:numCombinations
                % Find right-most index to increase
                % j = find(ind < n-k+1:n, 1, 'last');
                for j = subsetSize:-1:1
                    if ind(j)<n-subsetSize+j
                        break;
                    end
                end

                % Increase index j, initialize all indices to j's right.
                % ind(j:k) = (ind(j) + 1) : (ind(j) + 1 + k - j);
                % P(i, :) = v(ind);
                for t=1:j-1
                    P(i, t) = predictoridx(ind(t));
                end
                indj = ind(j) - j + 1;
                for t = j:subsetSize
                    ind(t) = indj + t;
                    P(i, t) = predictoridx(indj + t);
                end
            end
        end

        function [modifiedModel, status] = convertModelToPredictOnMatrix(mdl)
            % table2MatrixForShapley converts a heterogeneous (containing multiple datatypes) table to a numeric matrix
            % For classreg models, making this change can result in a HUGE performance improvement
            % For example, for a classification tree trained on the adult dataset with 100 samples, and 14 predictors, a 20X speed improvement was observed
            % For ensembles this improvement was close to 10X.
            % This is probably because indexing operations are much faster for matrices
            % Further, the type inference need not be done for matrices
            % However, since the models were trained on heterogeneous tables, these models cannot directly predict on numeric matrices
            % To enable prediction on numeric matrices, these models must be reconstructed with the VariableRange in their DataSummary modified
            % Calling constructors is cumbersome because (1) classreg models do not have the same signature for constructors, and
            % (2) some of the inputs passed to these constructors are not already available in the trained model that needs to be reconstructed
            % The best solution is to use the codegen utilities toStruct and fromStruct (which run in MATLAB), with categorical data as well as 'categorical' type labels supported for codegen in 21a
            dataSummary = mdl.DataSummary;
            try
                % an earlier version of this software used
                % table2FitMatrix, however it is recommended to use the
                % information already in the model instead, to account
                % for the edge case that the model has seen more
                % categories in training than the X that is passed by
                % the user, this is useful for models that dummify
                % based on the total number of categories present in
                % the dataSummary e.g. svm
                variableRangeCategorical = dataSummary.VariableRange;
                numpreds = numel(variableRangeCategorical);
                variableRangeNumeric = cell(size(variableRangeCategorical)); % though the orientation row or column should not matter, make sure it is stored as a row just like it is in classreg models
                for idx = 1:numpreds
                    vrangepred = variableRangeCategorical{idx};
                    numCategories = numel(vrangepred);
                    if numCategories == 0
                        variableRangeNumeric{idx} = vrangepred; % empty
                    else
                        variableRangeNumeric{idx} = (1:numCategories)'; % encode as the model would have
                    end
                end
                % for gam, use its constructor, for all other models, call to/fromStruct methods
                if isa(mdl, 'ClassificationGAM')
                    % ClassificationGAM for example does not have a
                    % toStruct/fromStruct, but it does allow post-fit
                    % parameters to be passed to reconstruct the
                    % object.
                    % Use this modified/retrained object that
                    % knows how to predict on numeric matrices for
                    % speed on table inputs classification kernel does
                    % not allow that, so cshapley uses tables for
                    % prediction using kernel classifier/regressor.
                    %
                    % Note that it is still faster to do this operation
                    % instead of using predict on a table in a loop.
                    % For 100 samples from the census1994 dataset for
                    % example, doing this operation is 25X faster than
                    % using a table inside.
                    %
                    % The constructor for the compact classes is
                    % protected, otherwise, using the compact
                    % constructor would make more sense
                    dataSummary.VariableRange = variableRangeNumeric;
                    modifiedModel = ClassificationGAM(mdl.PrivX, mdl.PrivY,mdl.W, mdl.ModelParams,...
                        dataSummary,mdl.ClassSummary,...
                        mdl.PrivScoreTransform);
                    modifiedModel = compact(modifiedModel);
                elseif isa(mdl, 'RegressionGAM')
                    dataSummary.VariableRange = variableRangeNumeric;
                    modifiedModel = RegressionGAM(mdl.PrivX, mdl.PrivY,mdl.W, mdl.ModelParams,...
                        dataSummary,mdl.PrivResponseTransform);
                    modifiedModel = compact(modifiedModel);
                elseif isa(mdl, 'ClassificationECOC')
                    dataSummary.VariableRange = variableRangeNumeric;
                    modifiedModel = ClassificationECOC(mdl.PrivX, mdl.PrivY,mdl.W, mdl.ModelParams,...
                        dataSummary,mdl.ClassSummary, mdl.PrivScoreTransform);
                    modifiedModel = compact(modifiedModel);
                else
                    % call toStruct to convert it to a struct
                    mdlst = mdl.toStruct;
                    % modify the dataSummary VariableRange
                    dataSummary = mdlst.DataSummary;
                    % do the same operations as classifToStruct to use a struct for
                    % Variable range because fromStruct expects a struct variable range
                    cellVariableRange = variableRangeNumeric;
                    for ii = 1:numel(cellVariableRange)
                        fn = strcat('X',num2str(ii));
                        dataSummary.VariableRange.(fn) = cellVariableRange{ii};
                    end
                    dataSummary.TableInput = false; % this model now predicts on numeric matrices
                    mdlst.DataSummary = dataSummary;
                    struct2object = str2func(mdlst.FromStructFcn);
                    modifiedModel = struct2object(mdlst); % convert it back to a classreg object
                end
                status = 1;
            catch % classification kernel does not support to/fromStruct, nor does it support the constructor call for a fitted object like say ClassificationSVM, ClassificationGAM
                status = -1;
                modifiedModel = mdl; % return the unmodified model
            end
        end

        function [istree, isens, hasSurrogateSplits, isCutVariable] = isCandidateForTreeSHAP(mdl, M)
            isclassiftree = @(mdl) isa(mdl,'classreg.learning.classif.CompactClassificationTree');
            isclassifens = @(mdl)  isa(mdl,'classreg.learning.classif.CompactClassificationEnsemble') &&...
                ~isempty(mdl.Trained) && ... % if non empty inspect the learners for trees
                (all(cellfun(@(x)isa(x, 'classreg.learning.classif.CompactClassificationTree'), mdl.Trained)) || all(cellfun(@(x)isa(x, 'classreg.learning.classif.CompactClassifByBinaryRegr'), mdl.Trained))) &&... % cshapley values exhibit linearity, non-linear transformations won't work
                mdl.Trained{1}.ScoreTransform=="none";

            isregrtree = @(mdl) isa(mdl,'classreg.learning.regr.CompactRegressionTree');
            isregrens = @(mdl)  isa(mdl,'classreg.learning.regr.CompactRegressionEnsemble') &&...
                ~isempty(mdl.Trained) && ... % if non empty inspect the learners for trees
                all(cellfun(@(x)isa(x, 'classreg.learning.regr.CompactRegressionTree'), mdl.Trained)) && ...
                mdl.Trained{1}.ResponseTransform=="none"; % cshapley values exhibit linearity, non-linear transformations won't work

            istree = isclassiftree(mdl) || isregrtree(mdl);
            isens = isclassifens(mdl) || isregrens(mdl);
            % look for surrogate splits and have a boolean for tracking
            % which variables are cutvar because if those columns in the
            % data have nans, we need to use interventional kernel
            isCutVariable = []; % default empty
            hasSurrogateSplits = false; % default false
            if istree
                hasSurrogateSplits = ~isempty(mdl.Impl.SurrCutVar);
            elseif isens
                hasSurrogateSplits = any(cellfun(@(x) ...
                    (isa(x, 'classreg.learning.classif.CompactClassificationTree') && ~isempty(x.Impl.SurrCutVar)) || ...
                    (isa(x, 'classreg.learning.regr.CompactRegressionTree') && ~isempty(x.Impl.SurrCutVar)) || ...
                    (isa(x, 'classreg.learning.classif.CompactClassifByBinaryRegr') && ...
                    ~isempty(x.CompactRegressionLearner.Impl.SurrCutVar)), mdl.Trained));
            end
            if hasSurrogateSplits
                isCutVariable = false(1,M);
                if istree
                    impl = treeimpl(mdl);
                    idx = unique(impl.CutVar);
                    idx = idx(idx>0);
                    isCutVariable(idx) = true;
                else
                    trees = mdl.Trained;
                    numTrained = numel(trees);
                    for i = 1:numTrained
                        t = trees{i};
                        impl = treeimpl(t);
                        idx = unique(impl.CutVar);
                        idx = idx(idx>0);
                        isCutVariable(idx) = true;
                    end
                end
            end
        end


        function islinear = isCandidateForLinearSHAP(mdl)
            hasNoOrdinalPredictors = @(mdl) ~isa(mdl, 'function_handle') && ~any(mdl.DataSummary.OrdinalPredictors); % linear shap cannot handle ordinal data
            isregrsvm = @(mdl) (isa(mdl, 'RegressionSVM') || isa(mdl,'classreg.learning.regr.CompactRegressionSVM')) && mdl.KernelParameters.Function=="linear";
            isclassifsvm = @(mdl) (isa(mdl, 'ClassificationSVM') || isa(mdl,'classreg.learning.classif.CompactClassificationSVM')) && mdl.KernelParameters.Function=="linear";
            isclassif = @(mdl) isa(mdl, 'ClassificationLinear') || isclassifsvm(mdl);
            isregr =  @(mdl) isa(mdl, 'RegressionLinear') || isregrsvm(mdl);

            islinear = isclassif(mdl) || isregr(mdl);
            islinearsvm = isregrsvm(mdl) || isclassifsvm(mdl);
            islinear = hasNoOrdinalPredictors(mdl) && (islinear || islinearsvm);
        end


        function tf = isTreeEnsembleGAM(mdl)
            % status is zero, so just need to check for compact models
            isclassif = @(mdl) isa(mdl,'classreg.learning.classif.CompactClassificationTree') || isa(mdl,'ClassificationTree') ... % tree
                || isa(mdl, 'ClassificationGAM') || isa(mdl,'classreg.learning.classif.CompactClassificationGAM')... % gam uses trees
                || isa(mdl,'classreg.learning.classif.CompactClassificationEnsemble') || isa(mdl,'classreg.learning.classif.ClassificationEnsemble') ... %
                && ~isempty(mdl.Trained) && ... % if non empty inspect the learners for trees
                (all(cellfun(@(x)isa(x, 'classreg.learning.classif.CompactClassificationTree'), mdl.Trained)) || all(cellfun(@(x)isa(x, 'classreg.learning.classif.CompactClassifByBinaryRegr'), mdl.Trained)));

            isregr = @(mdl) isa(mdl,'classreg.learning.regr.CompactRegressionTree') || isa(mdl,'RegressionTree')...
                || isa(mdl, 'classreg.learning.regr.CompactRegressionGAM') || isa(mdl, 'RegressionGAM')...
                || isa(mdl,'classreg.learning.regr.CompactRegressionEnsemble') || isa(mdl,'classreg.learning.regr.RegressionEnsemble')...
                && ~isempty(mdl.Trained) && ... % if non empty inspect the learners for trees
                all(cellfun(@(x)isa(x, 'classreg.learning.regr.CompactRegressionTree'), mdl.Trained));
            tf = isclassif(mdl) || isregr(mdl);
        end

        function checkForOrdinalsInData(thisX,isOrdinal)
            % this additional check is only performed for
            % conditional-kernel because distances cannot be computed for
            % such data
            % also note that models such as trees, gam do not identify
            % ordinal predictors in the data summary so the data needs to
            % be further probed for ordinals (the else branch)
            if isOrdinal
                error(message('stats:responsible:shapley:OrdinalsNotSupportedForConditional'));
            else % for models such as gam, trees OrdinalPredictors are not identified in DataSummary
                x = thisX(1,:); % thisX is non-empty for function handles, just take one row out for type inference
                if istable(x) % ordinals are only possible in tables, ensure non-empty args
                    [~,~,~,~,args] = classreg.learning.internal.table2FitMatrix(x,[]); % for a matrix, args is empty
                    if any(args{6}) % these are the ordinal predictor indices
                        error(message('stats:responsible:shapley:OrdinalsNotSupportedForConditional'));
                    end
                end
            end
        end

        function [classNameDefault,cnames] = defaultClassNameForShapleyPlot(mdl, isClassification)
            if isClassification
                cnames = mdl.ClassSummary.ClassNames;
                classNameDefault = cellstr(cnames(1));
            else
                classNameDefault = [];
                cnames = [];
            end
        end

        function [numpreds,classidx] = validateNumImportantPredictorsAndClassNameForSwarmBoxChart(isClassification,userProvided,numpreds,cname,cnames,M)
            validateattributes(numpreds, {'numeric'},...
                {'scalar','positive','finite', 'integer','<=',M});
            if ~isClassification
                if userProvided.ClassName % if the user specifies ClassNames
                    error(message('stats:responsible:shapley:ClassNamesInvalidForRegression'));
                end
                classidx = 1; % for regression and function handles there is only one column
            else % Classification
                userClassName = classreg.learning.internal.ClassLabel(cname);
                [~,classidx] = ismember(userClassName,cnames);
                if any(classidx==0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
                end
            end
        end

        function [predictorName, className, colorPredictorName, predictorIdx, classIdx, colorPredictorIdx] = validateArgsForDependencePlot( ...
                isClassification, userProvided, predictor, selectedClass, ...
                colorPredictor, predictorNames, classNames)

            validateattributes(predictor, {'string','char','numeric'}, {'nonempty'}, 1);
            if isnumeric(predictor)
                if predictor <= length(predictorNames) && ismember(predictor,1:length(predictorNames))
                    predictorIdx = predictor;
                else
                    predictorIdx = 0;
                end
            else
                [~, predictorIdx] = ismember(predictor, predictorNames);
            end
            if predictorIdx == 0
                error(message('stats:responsible:shapley:DependencePlotInvalidPredictorInput'));
            else
                predictorName = predictorNames(predictorIdx);
            end

            validateattributes(colorPredictor, {'string','char','numeric'},{},'','ColorPredictor');
            if isempty(colorPredictor)
                colorPredictorName = '';
                colorPredictorIdx = [];
            else
                if isnumeric(colorPredictor) 
                    if colorPredictor <= length(predictorNames) && ismember(colorPredictor, 1:length(predictorNames))
                        colorPredictorIdx = colorPredictor;
                    else
                        colorPredictorIdx = 0;
                    end
                else
                    [~, colorPredictorIdx] = ismember(colorPredictor, predictorNames);
                end
                if colorPredictorIdx == 0
                    error(message('stats:responsible:shapley:DependencePlotInvalidPredictorInput'));
                else
                    colorPredictorName = predictorNames(colorPredictorIdx);
                end
            end

            if ~isClassification
                if userProvided.ClassName % if the user specifies ClassNames
                    error(message('stats:responsible:shapley:ClassNamesInvalidForRegression'));
                end
                classIdx = 1; % for regression and function handles there is only one column
                className = '';
            else % Classification
                className = classreg.learning.internal.ClassLabel(selectedClass);
                [~,classIdx] = ismember(className, classNames);
                if any(classIdx==0)
                    error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
                end
            end
        end
    end
end

function [cutvar, children, cutpoint, isbranch, cutcats, response] = treeParameters(tree, isclassif)
treeimpl = tree.Impl;
if isclassif
    response = treeimpl.ClassProb;
else
    response = treeimpl.NodeMean;
end
response = response'; % transpose for built in
cutvar = double(treeimpl.CutVar); % always double
children = double(treeimpl.Children'); % always double and transpose for builtin
cutpoint = treeimpl.CutPoint;
isbranch = treeimpl.IsBranch;
cutcats = treeimpl.CutCategories;
end

function impl = treeimpl(t)
if isa(t, 'classreg.learning.classif.CompactClassifByBinaryRegr')
    t = t.CompactRegressionLearner;
end
impl = t.Impl;
end

function useKernelSHAP = checkForSurrogates(hasSurrogateSplits, Xmat, qmat,isCutVariable,isBatchFit)
if isBatchFit
    n = size(qmat,1);
    useKernelSHAP = false(n,1);
else
    useKernelSHAP = false;
end
if hasSurrogateSplits % do deeper checks to determine what methods to use for surrogate splits
    nanPredictorsInReference = any(isnan(Xmat),1);
    if any(isCutVariable & nanPredictorsInReference)
        % this means we cannot use tree shap
        useKernelSHAP = ~useKernelSHAP;
    else % reference predictors do not have any nan (and also in the tree)
        if isBatchFit
            nansInQuery = isnan(qmat);% qmat could still have nans
            useKernelSHAP = any((nansInQuery & isCutVariable),2); % use interventional tree for the query points you can
        else % one query point
            nanPredictorsInQuery = isnan(qmat);
            useKernelSHAP = any(nanPredictorsInQuery(:) & isCutVariable(:));
        end
    end
end
end

% Helper to determine whether a model is tree based, and if it is, whether
% the model does NOT contain objects (for example, gpuArrays).
function tf = iIsModelTreeBasedAndDoesNotContainObjects(mdl, istree, isens)
if istree
    treeToInspect = mdl;
elseif isens
    % Note that this assumes that the ensemble consists of learners of the
    % same type of tree.
    if isa(mdl.Trained{1}, "classreg.learning.classif.CompactClassifByBinaryRegr")
        treeToInspect = mdl.Trained{1}.CompactRegressionLearner;
    else
        treeToInspect = mdl.Trained{1};
    end
else
    % The model is not a tree or an ensemble of trees, just return
    % false.
    tf = false;
    return
end

% If the model was fitted using objects, then CutPoint will be an object.
tf = ~isobject(treeToInspect.CutPoint);
end

% Helper method to validate color map argument for Shapley plots
function colorMap = iValidateColorMapArg(cmap)
    validateattributes(cmap, {'string','char','numeric'},{'nonempty'});

    if ischar(cmap) || isstring(cmap)
        colorMap = validatestring(cmap, ["default","parula","turbo", ...
            "hsv","hot","cool","spring","summer","autumn","winter",...
            "gray","bone","copper","pink","sky","abyss","jet",...
            "lines", "colorcube","prism","flag","white", "bluered"]);
    else % use the same validation as core matlab
        if size(cmap,2) ~= 3
            error(message('MATLAB:colormap:InvalidNumberColumns'));
        end

        if ~isa(cmap,'uint8') && (min(cmap,[],'all') < 0 || max(cmap,[],'all') > 1)
            error(message('MATLAB:colormap:InvalidInputRange'))
        end
        colorMap = cmap;
    end
end

% LocalWords:  BLACKBOX QUERYPOINT Blackbox blackbox BLACKBOXMODEL BLACKBOXFITTED SHAPLEYVALUES
% LocalWords:  NUMSUBSETS SHAP Aas Func Preds BBand NVP SHapley Planations XLabel XSlow qmat Xmat
% LocalWords:  shap categcols heteropdist scalartext EXPECTEDVALUES numsubsets numobservations
% LocalWords:  heterpdist goodall falses Enu j's gam SHAPVALS BATCHFIT QUERYPOINTS batchfit
% LocalWords:  outputfcn cutpoints NShapley isens swarmchart boxchart MEANABSOLUTESHAPLEY
% LocalWords:  RAWSHAPLEYVALUES YJitter parula bluered numpreds yjitter cmap jitteroutliers
% LocalWords:  numimportantpredictors varfun
