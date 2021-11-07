function funcao(c2, c3, c4, teta1, teta2) {
	return [
		(2*(c3**2))+(c2**2)+(6*(c4**2))-1, 
		(8*(c3**3)) + (6*c3*(c2**2)) + (36*c3*c2*c4) + (108*c3*(c4**2)) - teta1,
		(60*(c3**4)) + (60*(c3**2)*(c2**2)) + (576*(c3**2)*c2*c4) + (2232*(c3**2)*(c4**2)) + (252*(c4**2)*(c2**2)) + (1296*(c4**3)*c2) + (3348*(c4**4)) + (24*(c2**3)*c4) + (3*c2) - teta2
		];
}

function jAssociada(c2, c3, c4) {
	let linha1 = [(2*c2), 4*c3, 12*c4];
	let linha2 = [(12*c3*c2) + (36*c3*c4), ((24*(c3**2)) + (6*(c2**2)) + (36*c2*c4) + (108*(c4**2)) ),((36*c2*c3)+(108*c3*2*c4))];
	let linha3 = [(60*2*c2*(c3**2))+(576*(c3**2)*c4)+(252*2*(c4**2)*c2)+(1296*(c4**3))+(24*c4*(3*(c2**2))) + 3,
		(60*(4*(c3**3))) + 2*60*c3*(c2**2) + 576*2*c3*c2*c4 + 2232*2*c3*(c4**2),
		(576*(c3**2)*c2) +2232*2*(c3**2)*c4 + 252*2*c4*(c2**2)+(1296*3*(c4**2)*c2) + 3348*4*(c4**3) + 24*(c2**3)
		];
	return [linha1, linha2, linha3];
}

function newton(teta1, teta2, iter, X0) {
	X0 = [1,0,0];
	let currentIter = 0;
	let J, F, tolk;
	while (currentIter <= iter) {
		J = jAssociada(X0[0],X0[1],X0[2]);
		F = funcao(X0[0],X0[1],X0[2], teta1, teta2);

		deltaX = resolverSistemaEDeterminanteLU(J, F, -1).x.map((item)=>(-item));
		X0 = somaVetores(X0, deltaX);
	
		tolk = norma(deltaX)/norma(X0);
		if (tolk <= (10**(-4))) {
			return {X0, currentIter};
		}
		currentIter++;
	}
	return {mensagem: "Método não convergiu", currentIter, X0}
}



function broyden(teta1, teta2, iter, X0) {
	X0 = [1,0,0];
	let currentIter = 0;
	let J, F, tolk, X1, B;
	while (currentIter <= iter) {
		if (currentIter === 0) {
			J = jAssociada(X0[0],X0[1],X0[2]);
		} else {
			X0 = X1;
			J = associadaB(J,Y,deltaX);;
		}
		F = funcao(X0[0],X0[1],X0[2], teta1, teta2);
		deltaX = resolverSistemaEDeterminanteLU(J, F, -1).x.map((item)=>(-item));
		X1 = somaVetores(X0, deltaX);
		Y = subtracaoVetores(funcao(X1[0],X1[1],X1[2], teta1, teta2),funcao(X0[0],X0[1],X0[2], teta1, teta2));
		//console.log({Y,X0,X1,deltaX})
		tolk = norma(deltaX)/norma(X1);
	
		if (tolk <= (10**(-4))) {
			return {X1, currentIter};
		}
		
		currentIter++;
	}
	return {mensagem: "Método não convergiu", X1, currentIter}
}

function associadaB(B,Y,deltaX) {
	deltaX = deltaX.map((item) => [item]);
	Y = Y.map((item) => [item]);
	//console.log({Y, deltaX})
	let deltaXTransposta = matrizTransposta(deltaX);
	let denominador = multiplicaMatrizes(deltaXTransposta, deltaX)[0][0];
	let numerador = multiplicaMatrizes(subtracaoMatrizes(Y, multiplicaMatrizes(B,deltaX)),deltaXTransposta);
	
	return somaMatrizes(B, multiplicaMatrizPorCte(numerador, (1/denominador)));
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

function somaMatrizes(A, B) {
	let matriz = [...A];
	for (let i = 0; i < A.length; i++){
		for (let j = 0; j < A[i].length; j++) {
			matriz[i][j] = A[i][j] + B[i][j];
		}
	}
	return matriz;
}

function subtracaoMatrizes(A, B) {
	let matriz = [...A];
	for (let i = 0; i < A.length; i++){
		for (let j = 0; j < A[i].length; j++) {
			matriz[i][j] = A[i][j] - B[i][j];
		}
	}
	return matriz;
}

function multiplicaMatrizPorCte(A, cte) {
	let matriz = [...A];
	for (let i = 0; i < A.length; i++){
		for (let j = 0; j < A[i].length; j++) {
			matriz[i][j] = A[i][j] * cte;
		}
	}
	return matriz;
}

function somaVetores(A, B) {
	let vetor = new Array(A.length);
	for (let i = 0; i < A.length; i++){
		vetor[i] = A[i] + B[i];
	}
	return vetor;
}

function subtracaoVetores(A, B) {
	let vetor = new Array(A.length);
	for (let i = 0; i < A.length; i++){
		vetor[i] = A[i] - B[i];
	}
	return vetor;
}

function norma(X) {
	let valorAcumulador = 0;
	for (let i = 0; i < X.length; i++){
		valorAcumulador += (X[i])**2
	}

	return Math.sqrt(valorAcumulador);
}

function alterasSinalDeMatriz(A) {
	let result = [...A];
	for (let i = 0; i < result.length; i++){
		for (let j = 0; j < result[i].length; j++){
			result[i][j] = A[i][j]*(-1)
		}
	}
	return result;
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
	let determinanteA = determinanteComLU(L, U);
	y = substituicaoParaFrente(L, b);
	x = retroSubstituicao(U, y);
	if (det > 0) {
		return {x, determinanteA};
	} else {
		return {x};
	}
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

console.log(newton(0.75,6.5,1000));
console.log(broyden(0.75,6.5,1000));