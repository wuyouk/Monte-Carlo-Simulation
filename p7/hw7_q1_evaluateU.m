function evalU = hw7_q1_evaluateU(x, y)
%
if (x > -1) && (x < 0) && (y > -1) && (y < 0)
    evalU = -log(1/2);
elseif (x > 0) && (x < 1) && (y > 0) && (y < 1)
    evalU = -log(1/2);
else
    evalU = 0;
end

