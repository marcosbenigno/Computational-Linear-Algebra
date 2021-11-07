function power(A, N, TOLm, IDET, tol, X0) {
	tol = tol ? tol : 10**(-3);
	X0 = X0 ? X0 : [];
	let mensagem = "O método não permite o cálculo do determinate de A.";
	let Y = [];
	let historicoTol= [];
	for (let i = 0; i < N; i++) {
		if (X0.length != N) {
			X0.push([1]);
		}
		Y.push([0]);
	}
	
	let iter = 0;
			let lambda0;
			let lambda1;
			let tolAtual;
	do {

		if (iter == 0) {
			lambda0 = X0[0];
		} else {
			lambda0 = lambda1;
			X0 = Y;
		}
	Y = multiplicaMatrizes(A, X0);
	
	lambda1 = Y[0][0];
	for (let i = 0; i < N; i++) {
		Y[i][0] = Y[i][0]/lambda1;
	}
	tolAtual = R(lambda0, lambda1);
	historicoTol.push(tolAtual);
	iter++;
	} while (tolAtual >= tol && iter <= TOLm);
	let lambda = lambda1;
	let X = Y;
	
	return {X, lambda, iter, historicoTol, mensagem};
}

function R(lambda0, lambda1) {
	return Math.abs(lambda1 - lambda0)/Math.abs(lambda1);
}

function multiplicaMatrizes(m1, m2) {
    var result = [];
    for (var i = 0; i < m1.length; i++) {
        result[i] = [];
        for (var j = 0; j < m2[0].length; j++) {
            var sum = 0;
            for (var k = 0; k < m1[0].length; k++) {
			//	console.log({k,j})
                sum += m1[i][k] * m2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

function matrizTransposta(m) {
	return m[0].map((x,i) => m.map(x => x[i]));
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

function jacobi(A, N, TOLm, IDET, tol) {
	tol = tol ? tol : 10**(-2);
	let mensagem = 'A matriz precisa ser simétrica para a convergência.\n';
	if (!eSimetrica) {
		mensagem += "A não é simétrica. O método de jacobi não fornecerá resultados corretos.\n";
	}
	mensagem += "A é simétrica.";
	let X = geraIdentidade(N);
	let A1;
	let P;
	let PTransposta;
	let iter = 0;
	let historicoTol = [];
	while(!elementosForaDaDiagonalMenorQue(A1 ? A1 : A, tol) && iter <= TOLm) {
		if (iter != 0){
			A = A1;
			
		}
		P = matrizAssociadaP(A);
		
		PTransposta = matrizTransposta(P);
		A1 = multiplicaMatrizes(multiplicaMatrizes(PTransposta, A), P);

		//console.log({P, A});
		X = multiplicaMatrizes(X, P);
		iter++;
	}
	
	if (IDET > 0) {
		let determinanteA = 1;
		for (let i = 0; i < X.length; i++) {
			
			determinanteA *= A1 ? A1[i][i] : A[i][i];
		}
		return { X, lambda: A1 ? A1 : A, iter, determinanteA, mensagem};
		
	}else {
	
		return { X, lambda: A1 ? A1 : A, iter, mensagem};
	}
}

function elementosForaDaDiagonalMenorQue(A, tol) {
	for (let i = 0; i < A.length; i++) {
		for (let j = 0; j < A[i].length; j++) {
			if ((i != j) && Math.abs(A[i][j]) > tol) {
				return false;
			}
		}
	}
	return true;
}

function matrizAssociadaP(A) {
	let { maiorValorAbsoluto, i, j } = maiorValorAbsolutoForaDaDiagonal(A);
	let P = geraIdentidade(A.length);
	let phi;
	if (A[i][i] != A[j][j]) {
		phi = 0.5*(Math.atan((2*A[i][j])/(A[i][i]-A[j][j])));
	} else {
		phi = Math.PI/4;
	}

	P[i][i] = Math.cos(phi);
	P[i][j] = -Math.sin(phi);
	P[j][i] = Math.sin(phi);
	P[j][j] = Math.cos(phi);
	return P;
}

function maiorValorAbsolutoForaDaDiagonal(A) {
	let temp;
	let iFinal, jFinal;
	for (let i = 0; i < A.length; i++) {
		
		for (let j = 0; j < A[i].length; j++) {
			if (i != j) {
				if (temp === undefined) {
					temp = Math.abs(A[i][j])
					iFinal = i;
					jFinal = j;
					
				} else if (Math.abs(A[i][j]) > temp) {
					temp = Math.abs(A[i][j]);
					iFinal = i;
					jFinal = j;
				}
			}
		}
	}
	return { maiorValorAbsoluto: temp, i: iFinal, j: jFinal };	
}

function geraIdentidade(N) {
	let indentidade = new Array();
	for (let i = 0; i < N; i++) {
		indentidade.push([]);
		for (let j = 0; j < N; j++) {
			if (i == j) {
				indentidade[i][j] = 1;
			} else {
				indentidade[i][j] = 0;
			}
		}
	}
	return indentidade;
}
