
function norma(X) {
	let valorAcumulador = 0;
	for (let i = 0; i < X.length; i++){
		valorAcumulador += (X[i])**2
	}

	return Math.sqrt(valorAcumulador);
}

function residuo(X0, X1) {
	let X1menX0 = new Array(X0.length);
	for (let i = 0; i < X0.length; i++){
		X1menX0[i] = X1[i] - X0[i];
	}

	return (norma(X1menX0))/(norma(X1));;
}

function eDiagonalDominante(A) {
	for (let i = 0; i < A.length; i++) {
		let acumuladorLinhas = 0;
		let acumuladorColunas = 0;
			for (let j = 0; j < A[i].length; j++) {
					if (i != j) {
						acumuladorLinhas += Math.abs(A[i][j]);
						acumuladorColunas += Math.abs(A[j][i]);
					}
			}
		if (acumuladorLinhas > Math.abs(A[i][i]) || acumuladorColunas > Math.abs(A[i][i])){
			return false;
		}
	}
	return true;
}






function eSimetrica(A) {
	for (let i = 0; i < A.length; i++) {
		for (let j = 0; j < A[i].length; j++) {
			if (A[i][j] != A[j][i]) {
				return false;
			}
		}
		
	}
	return true;
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


 function jacobi(A, b, tolM, X0, tol) {
	tol = tol ? tol : 10**(-3);
	X0 = X0 ? X0 : new Array(b.length).fill(1);
	let mensagem = '';
	if (eDiagonalDominante(A)) {
		mensagem += "A é diagonal dominante. Há garantia de convergência.\n";
	} else {
		mensagem += "A não é diagonal dominante. Não há garantia de convergência.\n";
	}
	mensagem += "Não é possível calcular a determinante com o método\n";
	let iter = 0;
	let X1 = new Array(b.length).fill(0);
	let historicoTol = [];
	let residuoAtual;
	do {
		if (iter > tolM) {
			break;
		}
		if (iter != 0) {
			X0 = X1.slice(0);
			X1 = X1.fill(0);

		}

		for (let i = 0; i < X1.length; i++) {
			X1[i] += b[i]/A[i][i];
				for (let j = 0; j < A[i].length; j++){
					if (i != j){
						X1[i] -= (A[i][j]*X0[j])/(A[i][i]);
					}
				}
		}
		iter += 1;
		residuoAtual = residuo(X0, X1);
		historicoTol.push(residuoAtual);
	} while (residuoAtual > tol);

	return {X1, iter, historicoTol, mensagem};
}


 function gaussSeidel(A, b,tolM, tol, X0) {
	tol = tol ? tol : 10**(-3);
	X0 = X0 ? X0 : new Array(b.length).fill(1);
	let mensagem = '';
	mensagem += "O método possui garantia de convergência se A for simétrica e definida positiva.\n";
	if (eSimetrica(A)){
		mensagem += "A é simétrica.\n";
	} else {
		mensagem += "A não é simétrica. Não há garantia de convergência.\n";
	} 
	mensagem += "Não é possível calcular a determinante com o método\n";
	console.log(A)
	let iter = 0;
	let X1 = new Array(b.length).fill(0);
	//console.log(X1)
	let historicoTol = [];
	let residuoAtual;
	do {
		if (iter > tolM) {
			break;
		}
		if (iter != 0) {
			X0 = X1.slice(0);
			X1 = X1.fill(0);

		}
		for (let i = 0; i < X1.length; i++) {
			X1[i] += b[i]/A[i][i];
				for (let j = 0; j < A[i].length; j++){
					if (i != j){
						if (j > i) {
							X1[i] -= (A[i][j]*X0[j])/(A[i][i]);
						} else {
							X1[i] -= (A[i][j]*X1[j])/(A[i][i]);
						}
					}
				}
		}
		iter += 1;
		residuoAtual = residuo(X0, X1);
		console.log({X0})
		historicoTol.push(residuoAtual);
	} while (residuoAtual > tol);
	
	return {X1, iter, historicoTol, mensagem};
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
	return {determinante: detL * detU, determinanteL: detL, determinanteU: detU};
}

 function resolverSistemaEDeterminanteLU(A, b, det) {
	let { L, U } = decomposicaoLU(A);
	let determinante = determinanteComLU(L, U);
	let determinanteA = determinante.determinante;
	y = substituicaoParaFrente(L, b);
	x = retroSubstituicao(U, y);
	if (det > 0) {
		return {x, determinanteA, determinanteL: determinante.determinanteL, determinanteU: determinante.determinanteU};
	} else {
		return {x};
	}
}

function decomposicaoCholesky(A) {
	let mensagem = '';
	mensagem += "O método possui bom funcionamento se A for simétrica e definida positiva.\n";
	if (eSimetrica(A)){
		mensagem += "A é simétrica.\n";
	} else {
		mensagem += "A não é simétrica. O método não possui bom funcionamento.\n";
	} 
	
	    let L = new Array(A.length).fill(0).map(
           x => new Array(A.length).fill(0));
		   
	for (let i = 0; i < A.length; i++) {
		L[i][i] = A[i][i];
		for (let k = 0; k < i; k++){
			L[i][i] -= (L[i][k])**2
		
		}
		L[i][i] = Math.sqrt(L[i][i]);
		for (let j = i+1; j < A.length; j++){

				L[j][i] = A[i][j]/L[i][i];
				
		for (let k = 0; k < i; k++){
			L[j][i] -= (L[i][k]*L[j][k])/L[i][i];

		
		}

				
		}
		
	}
	let U = matrizTransposta(L);
	return {L, U, mensagem};
	
}

 function resolverSistemaEDeterminanteCholesky(A, b, det) {
	let { L, U, mensagem } = decomposicaoCholesky(A);
	let determinante = determinanteComLU(L, U);
	let determinanteA = determinante.determinante;
	y = substituicaoParaFrente(L, b);
	x = retroSubstituicao(U, y);
	if (det > 0) {
		return {x, determinanteA, determinanteL: determinante.determinanteL, determinanteLTransposta: determinante.determinanteU, mensagem};
	} else {
		return {x, mensagem};
	}
}

function matrizTransposta(m) {
	return m[0].map((x,i) => m.map(x => x[i]));
}


