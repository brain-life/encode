% This function return the indices to the atoms having a particular spatial orientation determined by a main_orient +- offest
%  Copyright (2015), Franco Pestilli (Indiana Univ.) - Cesar F. Caiafa (CONICET)
%  email: pestillifranco@gmail.com and ccaiafa@gmail.com
function [ ind] = feGetAtoms(fe, main_orient, offset)
% INPUTS:
% fe: fe structure
% main_orient: [x,y,z] main orientation vector having unit-norm
% offset: maximum degrees appart from main orientation

% OUTPUT:
% ind: indices to atoms in the Dictionary meeting the criterion

ang = 180*acos(abs(main_orient*fe.life.M.orient))/pi;

ind = find(ang < offset);

end

