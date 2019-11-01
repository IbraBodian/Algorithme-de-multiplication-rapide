fft(ni,f,w)={\\ n (une puissance de 2) est la taille du tableau f avec f[1]=coeff cst et f[n]=coeff de n-1 et w est une prim root n-ème de 								\\l'unité
	n=#f;
	if(n==1,return([f[1]]));
	n=n/2;
	my(R0=vector(n));
	my(R1=vector(n));
	for( j=1,#R0,R0[j]=f[j]+f[j+n];R1[j]=(f[j]-f[j+n])*w^(j-1) );
	my( resultat=concat(fft(ni,R0,w^2),fft(ni,R1,w^2)) ) ;
	n=2*n;
	if(n!=ni,return(resultat),
		for(j=0,ni-1,
			L=Vecrev(digits(j,2));
			for(i=1,valuation(ni,2)-#L,L=concat(L,0));
			i=evaluation(Vecrev(L),2); \\fromdigits(L,2)
			if(i>=j,
				tmp=resultat[j+1];
				resultat[j+1]=resultat[i+1];
				resultat[i+1]=tmp;
			);
		);
		return(resultat);
	);
}

\\print(fft([1,1,1,1,1,1,1,0],Mod(4,257))); \\le polynome correspond à Polrev du vecteur en paramètre

evaluation(c,d)={
	return(sum(i=1,#c,c[i]*d^(i-1)));
}

fastconv(f,g,w)={
	n=#f;
	alpha=fft(n,f,w);
	beta=fft(n,g,w);
	gam=vector(#f);
	for( i=1,#f,gam[i]=alpha[i]*beta[i] );
	return(fft(n,gam,w^-1)/n);
}

\\print(fastconv([1,1,1,1,1,1,1,0],[1,1,1,1,1,1,0,1],Mod(4,257)));

findroot(k,m,n)={ \\p de la forme k*2^m+1 , on veut une racine primitive 2^n-eme
	gen=znprimroot(k*2^m+1);
	return(lift(gen^(k*2^(m-n))));
}

threeprime(a,b)={

	P=[71*2^57+1,75*2^57+1,95*2^57+1];
	W=[Mod(1830715477503673905,P[1]),Mod(545411398142084318,P[2]),Mod(12094553878055189361,P[3])]; \\racines 2^17-emes de l'unité
	A=Vecrev(digits(a,2^64));\\transforme a en poly
	B=Vecrev(digits(b,2^64));\\transforme b en poly
	ta=#A;
	tb=#B;
	t=ceil(log(ta+tb-1)/log(2));
	for(i=1,(2^t)-ta,A=concat(A,0));
	for(i=1,(2^t)-tb,B=concat(B,0));
	C=vector(3);
	for(j=1,3,
		AA=Mod(A,P[j]);
		BB=Mod(B,P[j]);
		C[j]=fastconv( AA,BB,W[j]^(2^(17-t)) );
	);
	CC=vector(2^t);
	for(j=1,2^t,CC[j]=lift( chinese( ( chinese(C[1][j],C[2][j]) ), C[3][j] ) ) );
	return(evaluation(CC,2^64));	
}

a=1684877070693736302020101041089974814184984649044849484968048746808466408068406468600406402849048903640084861000130337938337910971351507346103407343017310390981933987939143097193087991341439079134079308810437910490030710603781791089701398790137817817801937809139084098405687680749874098373893460749330980746848306364039689303008048634048930058684684908046989106746101727035815769480679637423132328340746488468468468346860024022060762072616601101031317331317137307103675169016941969713071606431003971376095461079164325301477674690146;

b=421219267673434075505025260272493703546246162261212371242012186702116602017101617150101600712262225910021215250032584484584477742837876836525851835754327597745483496984785774298271997835359769783519827202609477622507677650945447772425349697534454454450484452284771024601421920187468524593473365187332745186712076591009922325752012158512232514671171227011747276686525431758953942370169909355783082085186622117117117086715006005515190518154150275257829332829284326775918792254235492428267901607750992844023865269791081325369418672557;

\\print(threeprime(a^10,b^10));
\\print(a^10*b^10);








