function [Signal] = Eval_AP_Phantoms(Phantom,AcqPars,t_eval)
k_fun = length(Phantom.Function_Set);
maxcon = zeros(1,k_fun);

p = inputParser;

addParameter(p,'Scale_Img',1,@isnumeric); %scales simulated voxels to same magnitude as measured ones
addParameter(p,'Delay_Time',0,@isnumeric); %measured/simulated bolus arrival time map (in seconds)
addParameter(p,'Function_Parameters',{'default',[]},@iscell); %optional parameter arguments for one or more of the phantom functions: e.g., uptake parameters for EMM

parse(p,Phantom.Parameter_Set{:});

Scale_Img = p.Results.Scale_Img;
T0 = p.Results.Delay_Time;
FnPars = p.Results.Function_Parameters;

S_valid = Scale_Img > 0;
T_valid = (t_eval > T0) & (abs(T0) > 0);

Dynamic_pts = S_valid & T_valid;

has_interior_parameters = ~isequal(FnPars{1},'default');
if has_interior_parameters
    n = length(FnPars);
    FnPars_label = cell(1,n);
    for i=1:2:n
        FnPars_label{i} = FnPars{i};
    end
end

for label = 0:k_fun
    if label == 0
        Signal = Phantom.Const_Image;
    else
        eval_set = (Phantom.Function_Labels == label) & (Dynamic_pts);

        ison = ~isempty(find(eval_set,1));
        if ison
            [~,maxcon(label)] = Phantom.Function_Set{label}([]);
            maxsig = con2sig(maxcon);

            % function parameter sets structured as multi-cell 
            % object with a name-value pair structure; FnPars is
            % set of interior function parameters with nested
            % name-value pair structure.
            if has_interior_parameters
                for i=2:2:n
                    d = length(FnPars{n});
                    FnPars_label{i} = FnPars{i};
                    s = size(eval_set);
                    for q=1:d
                        sq = size(FnPars{i}{q});
                        if isequal(s,sq)
                            FnPars_label{i}{q} = FnPars{i}{q}(eval_set);
                        end
                    end
                end
            end

            Eval_Times = t_eval - (T0(eval_set));
            if has_interior_parameters
                Concentration = Phantom.Function_Set{label}(Eval_Times,FnPars_label{:});
            else
                Concentration = Phantom.Function_Set{label}(Eval_Times);
            end
            Concentration = Concentration.*Dynamic_pts(eval_set);

            Dynamic_Raw = con2sig(Concentration,...
                'FA',AcqPars.FlipAngle,'TR',AcqPars.TR...
                ,'TE',AcqPars.TE);
            Dynamic = Scale_Img(eval_set).*Dynamic_Raw/maxsig(label);
            Signal(eval_set) = Signal(eval_set) + Dynamic;
        end
    end
end

end