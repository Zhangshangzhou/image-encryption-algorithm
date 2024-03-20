function mtr = izigzagScan(vec, m, n)

% izigzagScan(vec, m, n) Inverse zigzag scan.

% mtr = izigzagScan(vec, m, n) inverse zigzag scan of vec, and the target

% matrix is m-by-n. The input parameter vec is an 1-by-k vector. If the

% length of vector is smaller than the size of mtr (i.e. k < m*n), then, the

% pargram will use 0 to fill the matrix. If k > m*n, the vector will be truncated.

%

% version: v1.0

% author: tianlan

% time: Dec 16, 2016

% CHECK INPUT ARGUMENTS

if nargin < 3

error('Not enough input arguments!');

elseif nargin > 3

error('Too many input arguments!');

end

% Get length of vector.

[mv,nv,tv] = size(vec);

if tv ~= 1 || mv ~=1

error('The input vector needs to be 1-by-k vector!');

end

if nv < m*n % if nv < m*n, then, use 0 to fill the matrix.

vec = [vec, zeros(1,m*n-nv)];

end

% REVERSE SCAN A VECTOR BY USING ZIGZAG ORDER.

mtr = zeros(m,n);

% zigzag scan.

x = 1;

y = 1;

vecNum = 0;

vecNum = vecNum + 1;

mtr(x,y) = vec(1,vecNum);

while 1

% arrive right-up border, needs to change direction.

if (y+1) <= n && x == 1 % towards right one step.

y = y + 1;

elseif (y+1) > n && (x+1) <= m % towards down one step.

x = x + 1;

else

break;

end

vecNum = vecNum + 1;

mtr(x,y) = vec(1,vecNum) ;

% judge the scan process is arrived the right-down corner (i.e. mtr(m,n)) or not.

if x == m && y == n

break;

end

% scan the matrix towards with left-down direction.

while (y-1) >= 1 && (x+1) <= m

y = y - 1;

x = x + 1;

vecNum = vecNum + 1;

mtr(x,y) = vec(1,vecNum);

end

% arrive the left-down border, needs to change direction.

if (x+1) <= m && y == 1% towards down one step.

x = x + 1;

elseif (y+1) <= n && x == m % towards right one step.

y = y + 1;

else

break;

end

vecNum = vecNum + 1;

mtr(x,y) = vec(1,vecNum);

% judge the scan process is arrived the right-down corner (i.e. mtr(m,n)) or not.

if x == m && y == n

break;

end

% scan the matrix towards with right-up direction.

while (x-1) >= 1 && (y+1) <= n

x = x - 1;

y = y + 1;

vecNum = vecNum + 1;

mtr(x,y) = vec(1,vecNum);

end

end