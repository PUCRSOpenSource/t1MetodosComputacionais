format long e
function newton( f, df, x0, tol, nmax)
	phi = (1 + sqrt(5))/2
	f = inline(f)
	df = inline(df)
	x(1,1) = 1
	x(2,1) = double(x0 - (f(x0)/df(x0)))
	x(3,1) = phi
	x(4,1) = abs(x(1)-x0)
	i = 2
	while xor(i <= nmax, tol > x(4,i-1))
		 x(1,i) = i
		 x(2,i) = x(2,i-1) - f(x(2,i-1))/df(x(2,i-1));
		 %x(2,i) = aux
		 x(3,i) = phi
		 x(4,i) = abs(x(2,i)-x(2,i-1))
		 i++
	end

	file = fopen('./tables/phi-newton.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & PHI & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência do número de ouro pelo método de Newton}\n')
	fprintf(file, '\\label{table:phi-newton}\n')
	fprintf(file, '\\end{table}')
	fclose(file)

end

function phi_frac(iteration=10, err=1*10^-1)
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
	fprintf(file,'Iteração & Aproximação & PHI & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência do número de ouro pelo método de frações continuadas}\n')
	fprintf(file, '\\label{table:phi-frac}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

function pi_sin(iteration=1, err=1*10^-1)
	
	x(1,1) = 1
	x(2,1) = 3 + sin(3)
	x(3,1) = pi
	x(4,1) = abs(pi - x(2,1))

	i = 2
	while xor(i <= iteration,  err > x(4, i-1))
		x(1,i) = i
		x(2,i) = x(i-1) + sin(x(i-1))
		x(3,i) = pi
		x(4,i) = abs(pi - x(i-1))
		i++
	end

	file = fopen('./tables/pi-sin.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & pi & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência de pi utilizando funções trigonométricas}\n')
	fprintf(file, '\\label{table:pi-sin}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

function pi_leibniz(iteration, err)
	x(1,1) = 1
	x(2,1) = 1
	x(3,1) = pi
	x(4,1) = abs(pi - x(2,1))

	i=2
	while xor(i <= iteration,  err > x(4, i-1))
		x(1,i) = i
		x(2,i) = x(i-1) + 2/((4*(i-1)-1)*4(i-1)+3)
		x(3,i) = pi
		x(4,i) = abs(pi - x(i-1))
		i++
	end
	file = fopen('./tables/pi-pow.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & pi & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência de pi utilizando funções trigonométricas}\n')
	fprintf(file, '\\label{table:pi-pow}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

function euler_taylor(iteration, err)
	x(1,1) = 1
	x(2,1) = 1
	x(3,1) = e
	x(4,1) = abs(e - x(2,1))
	i=2
	while xor(i <= iteration,  err > x(4, i-1))
		x(1,i) = i
		x(2,i) = x(i-1) + 1/factorial(i-1)
		x(3,i) = e
		x(4,i) = abs(e - x(i-1))
		i++
	end
	file = fopen('./tables/euler-taylor.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & e & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência de euler utilizando série de Taylor}\n')
	fprintf(file, '\\label{table:euler-taylor}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

function erdos(iteration, err)
	x(1,1) = 1
	x(2,1) = 1
	x(3,1) = abs(1.60669515241529176378330152319092458048057967150575643577807955369 - x(2,1))

	i=2
	while xor(i <= iteration,  err > x(3, i-1))
		x(1,i) = i
		x(2,i) = x(i-1) + 1/factorial(i-1)
		x(2,i) = x(i-1) + (1 / ((2^i)-1));
		x(3,i) = abs(x(2,i) - x(2,i-1))
		i++
	end
	file = fopen('./tables/erdos.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência de Erdős-Borwein}\n')
	fprintf(file, '\\label{table:erdos}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

function exponential(p, iteration, err)
	x(1,1) = 1
	x(2,1) = 1
	x(3,1) = e^p
	x(4,1) = abs(e^p - x(2,1))

	i=2
	while xor(i <= iteration,  err > x(4, i-1))
		x(1,i) = i
		x(2,i) = x(2,i-1) + (power(p,i-1)) / (factorial(i-1))
		x(3,i) = e^p
		x(4,i) = abs(x(3,i) - x(2,i-1))
		i++

	end
	file = fopen('./tables/exponential.tex', 'w')
	fprintf(file, '\\begin{table}[H]\n')
	fprintf(file, '\\centering \n')
	fprintf(file, '\\begin{tabular}{|c|c|c|c|}\n')
	fprintf(file, '\\hline \n')
	fprintf(file,'Iteração & Aproximação & $e^x$ & Erro \\\\ \n');
	fprintf(file, '\\hline \n')
	fprintf(file,'%d & %.14e &  %.14e & %.14e \\\\ \n\\hline\n',x);
	fprintf(file, '\\end{tabular}\n')
	fprintf(file, '\\caption{Convergência de $e^x$ utilizando série de Taylor}\n')
	fprintf(file, '\\label{table:erdos}\n')
	fprintf(file, '\\end{table}')
	fclose(file)
end

%newton('x^2 - x - 1', '2*x - 1', 30, 10^-15, 50)
%phi_frac(50, 10^-15)
%pi_sin(50, 10^-15)
%euler_taylor(50, 10^-15)
%exponential(5, 50, 10^-15)
%erdos(50, 10^-15)
