function funcao(t, y, yPDerivada, c, m, k, a1, a2, a3, w1, w2, w3) {
	return(-c * yPDerivada - k * y + a1 * Math.sin(w1 * t)+ a2 * Math.sin(w2 * t) + a3 * Math.sin(w3 * t)) * (1 / m);
}

function rungeKuttaNystrom(c, m, k, a1, a2, a3, w1, w2, w3,     h, tempoTotal) {
	
	let tk = 0;
	let y = 0;
	let yPrimeiraDerivada = 0;
	let N = Math.floor(tempoTotal/h);
	let allResults = [];
	let k1,k2,k3,k4,l,q;
	for (let i = 0; i < N; i++){
		tk=(i + 1) * h;
		k1 = (h / 2) * funcao(tk, y, yPrimeiraDerivada, c, m, k, a1, a2, a3, w1, w2, w3);
		q = (h / 2)* (yPrimeiraDerivada + k1 / 2);
		k2 = (h / 2)* funcao(tk + h / 2, y + q, yPrimeiraDerivada + k1, c, m, k, a1, a2, a3, w1, w2, w3);
		k3 = (h / 2)* funcao(tk + h / 2, y + q, yPrimeiraDerivada + k2, c, m, k, a1, a2, a3, w1, w2, w3);
		l = h * (yPrimeiraDerivada + k3);
		k4 = (h/ 2)* funcao(tk + h, y + l, yPrimeiraDerivada + 2*k3, c, m, k, a1, a2, a3, w1, w2, w3);
		allResults.push({t: tk, y, yPrimeiraDerivada, ySegundaDerivada: funcao(tk, y, yPrimeiraDerivada, c, m, k, a1, a2, a3, w1, w2, w3)});
		if (i != N-1) {
			y += h*(yPrimeiraDerivada + (k1 + k2 + k3) / 3);
			yPrimeiraDerivada += (k1 + 2 * k2 + 2 * k3 + k4) / 3;
		}
	}
	return {t: tk, y, yPrimeiraDerivada, ySegundaDerivada: funcao(tk, y, yPrimeiraDerivada, c, m, k, a1, a2, a3, w1, w2, w3), allResults};
	
}
console.log(rungeKuttaNystrom(1,1,1,1,0,0,1,1,1,0.1,100));