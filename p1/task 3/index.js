function regressaoLinear(X, Y, IDET, x) {
	let P = [];
	for (let i = 0; i < X.length; i++) {
		P.push([]);
		P[i][0] = 1;
		P[i][1] = X[i]
	}
	Y = Y.map((elm)=> [elm]);
	let PTransposta = matrizTransposta(P);
	let A = multiplicaMatrizes(PTransposta,P);
	let C = multiplicaMatrizes(PTransposta,Y);
	C = C.map((elm) => elm[0]);
	let resultados = resolverSistemaEDeterminanteLU(A, C, IDET);
	let yEstimado = x ? resultados.x[0] + resultados.x[1]*x : "x desejado não foi fornecido.";
	
	return { ...resultados, formattedExpression: `f(x) = ${resultados.x[0]} + ${resultados.x[1]}x`, yEstimado };
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

function resolverSistemaEDeterminanteLU(A, b, IDET) {
	let { L, U } = decomposicaoLU(A);
	y = substituicaoParaFrente(L, b);
	x = retroSubstituicao(U, y);
	if (IDET > 0) {
		let determinanteA = determinanteComLU(L, U);
		return {L, U, x, determinanteA};
	} else {
		return {L, U, x};
	}
}

function matrizTransposta(m) {
	return m[0].map((x,i) => m.map(x => x[i]));
}


function interpolacaoLagrange(X, Y, IDET, x) {
	let phis = new Array(X.length).fill(1);
	let y = 0;
	let mensagem = "Não é possível calcular determinantes com o método."
//	console.log(phis);
	for (let i = 0; i < X.length; i++) {
		for (let j = 0; j < X.length; j++) {
			if (j != i) {
				phis[i] *= (x-X[j])/(X[i]-X[j]);
			}
		}
		y += phis[i]*Y[i]
	}
	return { y, mensagem };
}
