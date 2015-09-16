output_precision(60)
function [x, ex] = newton( f, df, x0, tol, nmax)
	f = inline(f);
	df = inline(df);
	x(1) = double(x0 - (f(x0)/df(x0)));
	ex(1) = abs(x(1)-x0);
	k = 2;
	while k <= nmax && ex(k-1) > tol
		 x(k) = double(double(x(k-1)) - double((f(x(k-1))/df(x(k-1)))));
		 ex(k) = abs(x(k)-x(k-1));
		 k = k+1;
	end
end

function [x] =  phi_frac(precision=1)
	aux = 1;
	for i= 1:precision
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

function [ef, e_vec] = euller_taylor(iteration)
	ef(1) = 1;
	e_vec(1) = e;
	for i = 2:iteration
		ef(i) =ef(i-1) + 1/factorial(i-1);
		e_vec(i) = e;
	end
end

function [er] = erdos(iteration)
	er = 0
	for i=1:iteration
		er += 1 / ((2^i)-1);
	end
end

