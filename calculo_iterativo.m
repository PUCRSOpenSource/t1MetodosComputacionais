format long g
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

function phi_frac(iteration=1, err=1*10^-1)
	phi = (1 + sqrt(5))/2
	aux = 1
	x(1,1) = 1
	x(2,1) = aux
	x(3,1) = phi
	x(4,1) = abs(phi - aux)

	i = 2
	while xor(i <= iteration,  err > x(4, i-1))
		aux = double(1 + 1/aux)
		x(1,i) = i
		x(2,i) = aux
		x(3,i) = phi
		x(4,i) = abs(phi - aux)
		i++
	end

	file = fopen('./tables/phi.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação &  PHI & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14f &  %.14f & %.14f \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\label{table:phi-frac}\n')
	fprintf(file, '\\caption{Convergência do número de ouro pelo método de frações continuadas}\n')
	fprintf(file, '\\end{table}')
	fclose(file)

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
