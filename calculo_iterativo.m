output_precision(60)
function [x] = newton( f, df, x0, tol, nmax)
	f = inline(f);
	df = inline(df);
	x(1,1) = double(x0 - (f(x0)/df(x0)));
	x(2,1) = abs(x(1)-x0);
	k = 2;
	while k <= nmax && e(2,k-1) > tol
		 x(1,k) = double(double(x(1,k-1)) - double((f(x(1,k-1))/df(x(1,k-1)))));
		 ex(2,k) = abs(x(1,k)-x(1,k-1));
		 k = k+1;
	end
end

function [x] =  phi_frac(iteration=1)
	aux = 1;
	for i= 1:iteration
		aux = double(1 + 1/aux);
		x(i) = aux;
	end
end

function [pif, pi_vec] = pi_it(iteration)
	pif(1) = 3 + sin(3);
	pi_vec(1) = pi
	for i = 2:iteration
		pif(i) = pif(i-1) + sin(pif(i-1));
		pi_vec(i) = pi;
		aux = pif(i);
	end
end
function [pif, pi_vec] = pi_it2(iteration)
	pif(1) = 4;
	pi_vec(1) = pi
	for i = 2:iteration
		pif(i) = pif(i-1) + power(-1, mod(i+1,2)) * 4/(2*i - 1);
		pi_vec(i) = pi;
		aux = pif(i);
	end
end

function [ef, e_vec] = euler_taylor(iteration)
	ef(1) = 1;
	e_vec(1) = e;
	for i = 2:iteration
		ef(i) =ef(i-1) + 1/factorial(i-1);
		e_vec(i) = e;
	end
end

function [er] = erdos(iteration)
	er(1) = 1;
	for i=2:iteration
		er(i) = er(i-1) + (1 / ((2^i)-1));
	end
end

function [ex] = exponential(x, iteration)
	ex2 = e^x;
	ex(1,1) = 1;
	ex(2,1) = ex2;
	ex(3,1) = ex2 - ex(1,1);
	for k=2:iteration
		ex(1,k) = ex(1,k-1) + (power(x,k-1) / factorial(k-1));
		ex(2,k) = ex2;
		ex(3,k) = ex2 - ex(1,k);
	end
end
