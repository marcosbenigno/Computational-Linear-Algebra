
function funcao(c1, c2, c3, c4, x) {
	return c3*(x**c4)+(c1*Math.exp(c2*x));
}

function derivadaFuncao(c1,c2,c3,c4,x) {
	return (c1*c2*Math.exp(c2*x))+(c3*c4*(x**(c4-1)));
}

function bissecao(c1, c2, c3, c4, a, b, iter, tol) {
	tol = tol ? tol : 10**(-4);
	iter = iter ? iter : 1000;
	let xi;
	let valorF;
	let currentIter = 0;
	while(Math.abs(b-a) > tol) {
		if (currentIter == iter) {
			return {mensagem: "Não convergiu", x: xi, iter: currentIter};
		}
		xi = (a+b)/2;
		valorF = funcao(c1, c2, c3, c4, xi);
		if (valorF > 0) {
			b = xi;
		} else {
			a = xi;
		}
		currentIter++;
	}
	return {x: xi, currentIter};
}


function newton(c1, c2, c3, c4, a, b, iter, tol) {
	tol = tol ? tol : 10**(-4);
	iter = iter ? iter : 1000;
	let currentIter = 0;
	let x0 = (a+b)/2;
	let tolk, x;
	while (currentIter <= iter) {
		if (currentIter != 0) {
			x0 = x;
		}
		x = x0 - (funcao(c1, c2, c3, c4, x0)/derivadaFuncao(c1, c2, c3, c4, x0));

		tolk = Math.abs(x - x0);
		if (tolk <= tol) {
			return {x, currentIter};
		}
		currentIter++;
	}
	return {mensagem: "Não convergiu", x, currentIter};
}
//-- acima, item 1.
//-- abaixo, item 2.
function polinomial(c1, c2, c3, c4, a, b, N) {
	//ok funcionando
	let vandermonde = new Array(N);
	let vandermondeResult = new Array(N);
	let delta = (b-a)/(N-1);
	let resultado = 0;
	let pontos = [];
	let pesosW;
	if (N == 1){
		pontos = [(a + b)/2];
	}else{
		for (let i = 0; i < N; i++){
			pontos.push(a + (i * delta))
		}
	}
	
	for (let i = 0; i < N; i++) {
		vandermondeResult[i] = ((b**(i+1))-(a**(i+1)))/(i+1);
		vandermonde[i] = [];
		for (let j = 0; j < N; j++) {
			vandermonde[i][j] = ((pontos[j])**i);
		}
	}
	pesosW = resolverSistemaEDeterminanteLU(vandermonde, vandermondeResult, -1).x.map((item)=>-item);
	for (let i = 0; i < N; i++){
		resultado += pesosW[i] * funcao(c1,c2,c3,c4,pontos[i])
	}
	return {resultado: resultado};
}



function gaussLegendre(c1, c2, c3, c4, a, b, N) {
	//ok
	if (N < 2 || N > 10) {
		return {resultado: "", mensagem: "N não está no itervalo entre 2 e 10."};
	}
	const PONTOS = [null,null,
		[-0.5773502691896257, 0.5773502691896257], 
		[0.0000000000000000, -0.7745966692414834, 0.7745966692414834], 
		[-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526], 
		[0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640], 
		[0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521], 
		[0.0000000000000000, 0.4058451513773972, -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, -0.9491079123427585, 0.9491079123427585], 
		[-0.1834346424956498, 0.1834346424956498, -0.5255324099163290, 0.5255324099163290, -0.7966664774136267, 0.7966664774136267, -0.9602898564975363, 0.9602898564975363], 
		[0.0000000000000000, -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261, -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904], 
		[-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472, -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845, -0.9739065285171717, 0.9739065285171717]
		];
	const PESOSW = [null,null,
		[1.0000000000000000, 1.0000000000000000], 
		[0.8888888888888888, 0.5555555555555556, 0.5555555555555556], 
		[0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538], 
		[0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891], 
		[0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704], 
		[0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697], 
		[0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873, 0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763], 
		[0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354], 
		[0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881]
		];
    let L = b - a;
    let result = 0;
    let x = [];
    
    for (let i = 0; i < N; i++) {
      x.push((a + b + (PONTOS[N][i] * L))/2);
	  
	  result += funcao(c1, c2, c3, c4, x[i]) * PESOSW[N][i];
    }
    result = ((result * L)/2);
    return {resultado: result};
}

function substituicaoParaFrente(L, b) {
	let y = new Array(b.length).fill(0);
	y[0] = b[0]/L[0][0];
	for (let i = 1; i < y.length; i++) {
		y[i] = b[i]/L[i][i];
		for (let j = 0; j < i; j++) {
			 y[i] -= (L[i][j]*y[j])/L[i][i];
		}
	}
	return y;
}

function retroSubstituicao(U, y) {
	let x = new Array(y.length).fill(0);
	x[x.length-1] = y[x.length-1]/U[x.length-1][x.length-1]
    
    let s = 0;
	for (let i = x.length - 2; i > -1; i--){
        c = y[i]       
		for (let j = i+1; j < x.length; j++){
            c = c - U[i][j]*x[j]
		}	
        x[i] = c/U[i][i]
	}
   
	return x;
}


function decomposicaoLU(A) {
	let mat = A;
    let L = new Array(A.length).fill(0).map(
           x => new Array(A.length).fill(0));
    let U = new Array(A.length).fill(0).map(
           x => new Array(A.length).fill(0));

 
    for(let i = 0; i < A.length; i++) {
        //U
        for(let k = i; k < A.length; k++) {
            let acumulador = 0;
            for(let j = 0; j < i; j++) {

                acumulador += (L[i][j] * U[j][k]);
			
			}

            U[i][k] = A[i][k] - acumulador;
			
        }
 
        //L
        for(let k = i; k < A.length; k++) {
            if (i == k) {
                L[i][i] = 1;
			} else {
                let acumulador = 0;
                for(let j = 0; j < i; j++) {
                    acumulador += (L[k][j] * U[j][i]);
					
				}
				L[k][i] = (A[k][i] - acumulador)/U[i][i];
	            }
        }
    }
	return {L,U};
}

function determinanteComLU(L,U) {
	let detL = 1;
	let detU = 1;
	for (let i = 0; i < L.length; i++) {
		detL *= L[i][i];
		detU *= U[i][i];
	}
	return detL * detU;
}

 function resolverSistemaEDeterminanteLU(A, b, det) {
	let { L, U } = decomposicaoLU(A);
 //console.log({L,U})
	let determinanteA = determinanteComLU(L, U);
	y = substituicaoParaFrente(L, b);
	x = retroSubstituicao(U, y);
	if (det > 0) {
		return {x, determinanteA};
	} else {
		return {x};
	}
}

// item 2 acima
// item 3 abaixo

function diferencas(c1, c2, c3, c4, x, deltax, metodo="frente") {
	if (metodo == "frente") {
		return {resultado: (funcao(c1, c2, c3, c4, x + deltax) - funcao(c1, c2, c3, c4, x))/deltax};
	} else if (metodo == "tras") {
		return {resultado: (funcao(c1, c2, c3, c4, x)-funcao(c1, c2, c3, c4, x - deltax))/deltax}
	} else if (metodo == "central") {
		return {resultado: (funcao(c1, c2, c3, c4, x + deltax) - funcao(c1, c2, c3, c4, x - deltax))/(2*deltax)}
	} else {
		return {mensagem: "Método não permitido.", resultado: ""};
	}
}
// intem 3 acima
// item 4 abaixo

function extrapolacaoDeRichard(c1, c2, c3, c4, x, deltax1, deltax2, metodo="frente", p=1) {
	let d1, d2;
	let q = deltax1/deltax2;
	let result;
	if (metodo == "frente") {
		d1 = (funcao(c1, c2, c3, c4, x + deltax1) - funcao(c1, c2, c3, c4, x))/deltax1;
		d2 = (funcao(c1, c2, c3, c4, x + deltax2) - funcao(c1, c2, c3, c4, x))/deltax2;
	} else if (metodo == "tras") {
		d1 = (funcao(c1, c2, c3, c4, x)-funcao(c1, c2, c3, c4, x - deltax1))/deltax1;
		d2 = (funcao(c1, c2, c3, c4, x)-funcao(c1, c2, c3, c4, x - deltax2))/deltax2;
		
	} else if (metodo == "central") {
		d1 = (funcao(c1, c2, c3, c4, x + deltax1) - funcao(c1, c2, c3, c4, x - deltax1))/(2*deltax1);
		d2 = (funcao(c1, c2, c3, c4, x + deltax2) - funcao(c1, c2, c3, c4, x - deltax2))/(2*deltax2);
	} else {
		return {mensagem: "Método não permitido.", resultado: ""};
	}
	result = d1 + ((d1-d2)/((q**(-p))-1));
	return {resultado: result};
}

// item 4 acima