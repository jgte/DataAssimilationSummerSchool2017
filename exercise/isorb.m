function [out,mss] = isorb(orb)
% ISORB(ORB) returns true if the input is a valid orb as defined by the
% orb_constructor function. Otherwise it returns false.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>

% NOTE: NO ORB_* FUNCTIONS CAN BE CALLED FROM HERE. EVERY ORB_* FUNCTION
%       CALLS ISORB, AND THE PROGRAM IS STUCK IN AN INFINITE LOOP.

% Defaults
% NOTE: These are duplicated in the orb_set* function
LIST_REF_FRAMES = {'trf','crf','hill'};
LIST_ORB_TYPES = {'sph','cart','oe'};

% Allocate outputs
out = false(size(orb));
mss = cell(size(orb));

% Lopp all elements
for i = 1:numel(orb)
    % check if it is a structure
    if ~isstruct(orb(i))
        if nargout>1; mss{i}='Input orb is not a structure.'; end
        continue;
    % check if the fieldnames match the constructor
    elseif any(~isfield(orb(i),fieldnames(orb_constructor)))
        if nargout>1; mss{i}='Input orb does not have all required fields.'; end
        continue;
    else
        if nargout>1; mss{i}=[]; end
        out(i)=true;
    end

    %% Now check the optional fields
    % check the ReferenceDate
    if ~isempty(orb(i).ReferenceFrame) && ~isvector(orb(i).ReferenceDate) || numel(orb(i).ReferenceDate) ~= 6
        if nargout>1; mss{i}='Invalid ''ReferenceDate'' field.'; end
        out(i) = false, continue;
    end
    % check the ReferenceDate
    if ~isempty(orb(i).ReferenceFrame) && ~any(strcmp(orb(i).ReferenceFrame,LIST_REF_FRAMES))
        if nargout>1; mss{i}='Invalid ''ReferenceFrame'' field.'; end
        out(i) = false, continue;
    end
    % check the Units
    if ~iscell(orb(i).Units) || length(orb(i).Units)~=3
        if nargout>1; mss{i}='Invalid ''Units'' field.'; end
        out(i) = false, continue;
    end
    % check the Type
    if ~isempty(orb(i).Type) && ~any(strcmp(orb(i).Type,LIST_ORB_TYPES))
        if nargout>1; mss{i}='Invalid ''Type'' field.'; end
        out(i) = false, continue;
    end
    % Check the Description
    if ~ischar(orb(i).Comment)
        if nargout>1; mss{i}='Invalid ''Comment'' field.'; end
        out(i) = false, continue;
    end
end
