#
#                                                                                               
# Auxiliar functions for braid monodromies and groups.                                          
# Developed by:                                                                                 
#
#
# Enrique Artal Bartolo
#     ...Departamento de Matematicas-IUMA
# Universidad de Zaragoza
# artal@unizar.es
#
# Last update: 19-05-2015
#
#
from sage.groups.finitely_presented import wrap_FpGroup
#
def LibreNorm(a):
	r"""
	Return the cyclic reduction of `a` as member of a free group.
	
	INPUT:
	
	- ``a`` -- word of a free group (default: None)
	
	
	OUTPUT:

	The cyclic reduction of the word


	EXAMPLES: 

	This example illustrates a simple use of this function

	::

		sage: F=FreeGroup(3)
		sage: a=F([1,2,-1])
		sage: LibreNorm(a)
		x1


	"""
	F=a.parent()
	n=F.rank()
	L=a.Tietze()
	if len(L)<1:
		return (F(1))
	primero=(L[0]+L[-1]==0)
	while primero:
		L=L[1:-1]
		primero=(L[0]+L[-1]==0)
	return (F(L))


def LibreConj(a,b):
	r"""
	Checks if `a` and `b` are conjugate in their free group.
	
	INPUT:
	
	- ``a`` -- word of a free group (default: None)
	
	- ``b`` -- word of a free group (default: None)
	
	OUTPUT:

	``True`` if the words are conjugate, ``False`` if not.


	EXAMPLES:

	This example illustrates a simple use of this function

	::

		sage: F=FreeGroup(3)
		sage: a=F([2,1,-3,2])
		sage: b=F([-3,2,2,1])
		sage: LibreConj(a,b)
		True


	"""
	F=a.parent()
	n=F.rank()
	g=LibreNorm(a)
	h=LibreNorm(b)
	L=list(g.Tietze())
	M=list(h.Tietze())
	l=len(L)
	if l!=len(M):
		return (False)
	distintas=True
	j=0
	while distintas and j<l:
		distintas=(L!=M)
		L=L[1:]+[L[0]]
		j=j+1
	return (not distintas)

def revertirF(a):
	r"""
	Realizes the anti-automorphism $\\varphi:\\mathbb{F}_n\\to\\mathbb{F}_n$
	such that $x_i\\mapsto x_{n-1-i}$, $i=0,1\\dots,n-1$.
	
	INPUT:
	
	- ``a`` -- word of a free group (default: None)
	
	OUTPUT:

	Reverts the image by $\\varphi$ of `a`.


	EXAMPLES:

	This example illustrates a simple use of this function

	::

		sage: F=FreeGroup(3)
		sage: a=F([2,1,-3,2])
		sage: revertirF(a)
		x1*x2*x0^-1*x1


	"""
	F=a.parent()
	n=F.rank()
	L=a.Tietze()
	L1=(v.sign()*(n+1-v.abs()) for v in L)
	return (F(L1))

def revertirT(b):
	r"""
	Realizes the anti-automorphism $\\varphi:\\mathbb{B}_n\\to\\mathbb{B}_n$
	such that $x_i\\mapsto x_{n-2-i}$, $i=0,1\\dots,n-2$.
	
	INPUT:
	
	- ``b`` -- a braid (default: None)
	
	OUTPUT:

	Returns the image of `b` by $\\varphi$.


	EXAMPLES:

	This example illustrates a simple use of this function

	::

		sage: F=BraidGroup(3)
		sage: a=F([2,1,-2,1])
		sage: revertirF(a)
		s0*s1^-1*s0*s1


	"""
	B=b.parent()
	n=B.strands()
	L=b.Tietze()
	L1=(v.sign()*(n-v.abs()) for v in L)
	return (B(L1))

def LibreTrenza(libre,trenza):
	r"""
	It defines another left action of the braid group $\\mathbb{B}_n$ on the free group 
	$\\mathbb{F}_n$. Unlike the defined by ``x*s``, in this case, the image of $x_i$ and $s_i$
	equals $x_{i+1}$ and the image of $x_{i+1}$ and $s_i$ equals $x_{i+1} x_i x_{i+1}^{-1}$.
	
	INPUT:
	
	- ``libre`` -- an element of the free group of some rank `n`

	- ``trenza`` -- a braid in `n`of strands	
	OUTPUT:

	Returns this differently normalized braid action.


	EXAMPLES:

	This example illustrates a simple use of this function

	::

		sage: F=FreeGroup(3)
		sage: B=BraidGroup(3)
		sage: a=F([2,1,-3,2])
		sage: s=B([2])
		sage: a*s
		x1*x2*x1^-1*x0*x2*x1^-1
		sage: LibreTrenza(a,s)
		x2*x0*x2*x1^-1


	"""
	g=libre
	s=trenza
	return (revertirF(revertirF(g)*revertirT(s)))

def conjtrenza(lista,n):
	r"""
	It converts a list ``lista`` containing two lists of integers into a braid. Each sublist corresponds via Tietze 
	with a braid in ``n`` strands, say $t_1,t_2$. The result is the braid $t_1 t_2 t_1^{-1}$.
	
	INPUT:
	
	- ``lista`` -- a list of two sublists of non-zero integers coding two braids

	- ``n`` -- the number of strands
	
	OUTPUT:

	The left conjugation of the first braid to the second one.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: conjtrenza([[1,3,2],[2,1]],4) 
		s0*s2*s1^2*s0*s1^-1*s2^-1*s0^-1


	"""
	B=BraidGroup(n)
	u0,v0=[B(a) for a in lista]
	return (u0*v0/u0)

def relstrenzaconj(lista,n):
	r"""
	This function is to be used to compute a fundamental group from a braid monodromy. A list ``lista`` contains
	Tietze lists of integers to produce two braids $t_1,t_2$. The central braid is positive and connected and involves,
	$r$ strands from $i$ to $i+r-1$. The result is a list of words in the free group of the form 
	$(x_j^{-1} x_j^{t_2} )^{t_1^{-1}}$ for $j=i,\\dots,i+r-2$.
	
	INPUT:
	
	- ``lista`` -- a list of two sublists of non-zero integers coding two braids, the second one, positive and connected.

	- ``n`` -- the number of strands
	
	OUTPUT:

	A list of words in the free group of rank `n`.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: relstrenzaconj([[1,2],[3]],4)
		[x0^-1*x3]
		sage: relstrenza([[1,2],[3]],4)
		[x0^-1*x3, x1^-1*x3*x0^-1*x1*x0*x3^-1, x2^-1*x3*x0^-1*x2*x0*x3^-1, x0*x3^-1]


	"""
	B=BraidGroup(n)
	F=FreeGroup(n)
	a,b=lista
	t1=B(a)^-1
	t0=B(b)
	bv=[v.abs() for v in b]
	m=min(bv)
	M=max(bv)
	rel0=[v^-1*LibreTrenza(v,t0) for v in F.gens()[m-1:M]]
	return ([LibreTrenza(v,t1) for v in rel0])

def relstrenza(lista,n):
	r"""
	This function is to be used to compute a fundamental group from a braid monodromy. A list ``lista`` contains
	Tietze lists of integers to produce two braids $t_1,t_2$. It returns the list of $x_j^{-1} x_j^{t_2 t_1 t_2^{-1}}$.
	
	INPUT:
	
	- ``lista`` -- a list of two sublists of non-zero integers coding two braids.

	- ``n`` -- the number of strands
	
	OUTPUT:

	A list of `n` words in the free group of rank `n`.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: relstrenzaconj([[1,2],[3]],4)
		[x0^-1*x3]
		sage: relstrenza([[1,2],[3]],4)
		[x0^-1*x3, x1^-1*x3*x0^-1*x1*x0*x3^-1, x2^-1*x3*x0^-1*x2*x0*x3^-1, x0*x3^-1]
	"""
	F=FreeGroup(n)
	t0=conjtrenza(lista,n)
	rel0=[v^-1*LibreTrenza(v,t0) for v in F.gens()]
	return (rel0)

def invertirlista(lista):
	r"""
	Given a Tietze list ``lista``, it provides the Tietze list of the inverse.
	
	INPUT:
	
	- ``lista`` -- a list of non-zero integers coding a word in a free group
	
	
	OUTPUT:

	The Tietze list of the inverse of ``lista``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: invertirlista([2,-1,3])
		[-3, 1, -2]
	"""
	return ([-_ for _ in reversed(lista)])


def cambio_rel(lista,nuevo,elim):
	r"""
	Given two Tietze lists ``lista`` (representing a group word $w$), ``nuevo`` (representing another word $u$), 
	and a positive integer ``elim`` (representing a generator $x_i$), we replace each occurrence 
	of ``elim`` or its opposite in ``lista`` by the elements in ``nuevo`` (or its reversed opposite). The generator
	$x_i$ equals the word $w$, where a new generator $y_i$ (or its inverse) appears exactly once.
	
	INPUT:
	
	- ``lista`` -- a list of non-zero integers
	- ``nuevo`` -- a replacing list of non-zero integers
	- ``elim`` -- the integer to be replaced.
	
	
	OUTPUT:

	A list representing the Tietze list of the word after the replacement of the generator.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: cambio_rel([1,3,-5,2,6,5],[2,5,3,-6],5)
		[1, 3, 6, -3, -5, -2, 2, 6, 2, 5, 3, -6]
	"""
	pal=[]
	inv_nuevo=invertirlista(nuevo)
	for i in lista:
		if i.abs()!=elim:
			pal+=[i]
		elif i==elim:
			pal+=nuevo
		elif i==-elim:
			pal+=inv_nuevo
	return (pal)


def cambio(grupo,elim_nuevo,elim,lista=[]):
	r"""
	Given a group ``grupo``, a generator $x_i$ (represented by its index $i\\equiv$``elim``), the list
	``nuevo`` represents a word to eliminate $x_i$ in terms of a new generator which will have the same index.
	As an option one can add a list ``lista`` of Tietze words representing elements of the group.
	
	INPUT:
	
	- ``grupo`` -- a SAGE finitely presented group
	- ``elim_nuevo`` -- a replacing list of non-zero integers
	- ``elim`` -- the integer to be replaced.
	- ``lista`` -- an optional list (default: ``[]``) of words to rewrite.
	
	
	OUTPUT:

	The finitely presented group in the new generators, and the translation
	of the elements in ``lista``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(3)
		sage: rlk=[F.gen(0)^2,F([2,3])^2]
		sage: G=F/rlk
		sage: G1=cambio(G,[2,-3],2)[0]
		sage: G1.relations()
		(x0^2, x1^2)

	"""
	rels=[_.Tietze() for _ in grupo.relations()]
	m=len(grupo.generators())
	inv_elim_nuevo=invertirlista(elim_nuevo)
	rlk=[cambio_rel(_,elim_nuevo,elim) for _ in rels]
	listares=[cambio_rel(_,elim_nuevo,elim) for _ in lista]
	F0=FreeGroup(m)
	listares=[F0(_).Tietze() for _ in listares]
	g=F0/rlk
	P=g.gap().PresentationFpGroup()
	P.TzSearch()
	P.TzSearch()
	P.TzSearch()
	P.TzSearchEqual()
	return ([wrap_FpGroup(P.FpGroupPresentation()),listares])

def cambio_elim(lista,pal,elim):
	r"""
	Given two Tietze lists ``lista`` (representing a group word $w$), ``pal`` (representing another word $u$),
	where the generator $x_i$ (represented by the positive integer ``elim``), does not appear. 
	We replace each occurence of ``lista`` by the word $u$; we shift by $-1$ the indices
	greater that $i$ since we deal with one less generator.
	
	INPUT:
	
	- ``lista`` -- a list of non-zero integers
	- ``pal`` -- a replacing list of non-zero integers where $\\pm$``elim``does not appear.
	- ``elim`` -- the integer to be replaced.
	
	
	OUTPUT:

	A list with the replacement of ``elim`` and the shifting of indices.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: cambio_elim([2,1-3,4],[-2,1],3)
		[2, -2, 3]
	"""
	p=[]
	for i in lista:
		if i.abs()<elim:
			p.append(i)
		elif i.abs()>elim:
			p.append(i.sign()*(i.abs()-1 ))
		elif i==elim:
			p+=pal
		elif i==-elim:
			p+=palinv
	return (p)

def eliminar(grupo,generador,lista=[]):
	r"""
	Given a group ``grupo``, a generator $x_i$ (represented by its index $i\\equiv$ ``generador``), the lista
	``nuevo`` represents a word to eliminate $x_i$ in terms of a new generator which will have the same index.
	As an option one can add a list ``lista`` of Tietze words representing elements of the group. If the generator
	can be eliminated, it produces a new presentation with that generator erased; if not, the group remains unchanged.
	
	INPUT:
	
	- ``grupo`` -- a SAGE finitely presented group
	- ``generador`` -- the integer representing the generator to be eliminated.
	- ``lista`` -- an optional list (default: ``[]``) of words to rewrite.
	
	
	OUTPUT:

	The finitely presented group with, eventually, one less generator, and the translation
	of the elements in ``lista``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(3)
		sage: rlk=[F.gen(0)^2,F([1])*F([2,3])^2]
		sage: G=F/rlk
		sage: G1=eliminar(G,1)[0]
		sage: G1.relations()
		((x1^-1*x0^-1)^4,)

	"""
	rels=[_.Tietze() for _ in grupo.relations()]
	aux=[[i.abs() for i in _] for _ in rels]
	ind=[j for j in range(len(rels)) if aux[j].count(generador)==1 ]
	if ind==[]:
		return (grupo)
	cnt=[len(rels[j]) for j in ind]
	m=min(cnt)
	rel=rels[ind[cnt.index(m)]]
	if generador in rel:
		rel=invertirlista(rel)
	j=rel.index(-generador)
	pal0=rel[j+_sage_const_1 :]+rel[:j]
	pal=[]
	for i in pal0:
		if i.abs()<generador:
			pal.append(i)
		elif i.abs()>generador:
			pal.append(i.sign()*(i.abs()-1 ))
	palinv=invertirlista(pal)
	newrels=[cambio_elim(_,pal,generador) for _ in rels]
	listares=[cambio_elim(_,pal,generador) for _ in lista]
	n=len(grupo.generators())
	F=FreeGroup(n-1 )
	listares=[F(_).Tietze() for _ in listares]
	rlk=[F(_) for _ in newrels]
	g=F/rlk
	P=g.gap().PresentationFpGroup()
	P.TzSearch()
	P.TzSearch()
	P.TzSearchEqual()
	return ([wrap_FpGroup(P.FpGroupPresentation()),listares])

def CyclicComm(lista):
	r"""
	The elements of the list ``lista`` belong to a group.
	
	INPUT:
	
	- ``lista`` -- a list of elements of a group.
	
	
	OUTPUT:

	The list of reduced commutators of all (but one) elements of the list with the product of all of them.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(4)
		sage: L=F.gens()
		sage: CyclicComm(L)
		[x1*x0*x1*x2*x3*x1^-1*x3^-1*x2^-1*x1^-1*x0^-1,
		x2*x0*x1*x2*x3*x2^-1*x3^-1*x2^-1*x1^-1*x0^-1,
		x3*x0*x1*x2*x3^-1*x2^-1*x1^-1*x0^-1]

	"""
	n=len(lista)
	if n<2:
		return ([])
	pr=prod(lista)
	res=[]
	for a in lista[1:]:
		res.append(a*pr/a/pr)
	return (res)

def GtoK(G,elto,K):
	r"""
	Given an element ``elto``  a group ``G``, and a quotient ``K`` of it, write down its image in $K$. This function writes
	it as a word in the generators of $K$.
	
	INPUT:
	
	- ``G`` -- A finitely presented group of GAP.
	- ``elto`` -- An element of `K` as word in `G.
	- ``K`` -- A subgroup of `G`, element of GAP.
	
	
	
	OUTPUT:

	The element ``elto`` written in $K$.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(2)
		sage: rel=[F([1,2])^3]
		sage: G=F/rel
		sage: a=G([2])^2
		sage: g=G.gap()
		sage: u=a.gap()
		sage: K=g.FactorGroupFpGroupByRels([u])
		sage: K=g.FactorGroupFpGroupByRels([u^2])
		sage: GtoK(g,u,K)
		x1^2

	"""
	epiG=G.EpimorphismFromFreeGroup()
	epiK=K.EpimorphismFromFreeGroup()
	pal=epiG.PreImagesRepresentative(elto)
	num=pal.TietzeWordAbstractWord(epiG.Source().GeneratorsOfGroup())
	return (epiK.Image(num.AbstractWordTietzeWord(epiK.Source().GeneratorsOfGroup())))

def GtoTietze(G,elto):
	r"""
	Given an element ``elto`` of a GAP group ``G`` this function provides a Tietze list.
	
	INPUT:
	
	- ``G`` -- A finitely presented group of GAP.
	- ``elto`` -- An element of `G`.
	
	
	
	OUTPUT:

	A Tietze list.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(2); rel=[F([1,2])^3]; G=F/rel
		sage: a=G([2])^2; g=G.gap(); u=a.gap()
		sage: GtoTietze(g,u)
		[ 2, 2 ]

	"""
	epiG=G.EpimorphismFromFreeGroup()
	pal=epiG.PreImagesRepresentative(elto)
	num=pal.TietzeWordAbstractWord(epiG.Source().GeneratorsOfGroup())
	return (num.sage())

def abelianizar(tz,m):
	r"""
	Given an Tietze list ``tz`` from a free group of rank ``m``, it returns the abelianization vector.
	
	INPUT:
	
	- ``tz`` -- A Tietze list.
	- ``m`` -- An integer.
	
	
	
	OUTPUT:

	An integer vector such that the $i$-entry is the algebraic sum 
	of the number of $\\pm i$ entries.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: abelianizar([2,4,-1,2,1],4)
		(0, 2, 0, 1)

	"""
	v=vector(ZZ,m)
	L=FreeModule(ZZ,m)
	B=L.basis()
	for i in tz:
		signo=i.sign()
		absoluto=i.abs()
		v=v+signo*B[absoluto-1]
	return (v)

def MatrizAbel(grupo):
	r"""
	The matrix of the abelianization of the relations of a f.p. group.
	
	INPUT:
	
	- ``grupo`` -- A finitely presented group.
	
	
	
	OUTPUT:

	An integer valued matrix with as many rows as relations and as many columns
	as generators. Each row represents the abelianization of the relation.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(4)
		sage: MatrizAbel(B)
		[ 1 -1  0]
		[ 0  0  0]
		[ 0  1 -1]

	"""
	m=len(grupo.gens())
	TZ=map(lambda v:v.Tietze(),grupo.relations())
	A=[]
	for tz in TZ:
		A.append(abelianizar(tz,m))
	return (Matrix(A))

def CambioVarSmith(matriz,smith,tt):
	r"""
	The matrix ``matriz`` is an invertible matrix whose size $n$ is related to the number of generators
	of a group. The matrix ``smith`` is a diagonal $m\\times n$ where $m$ is related to the number of generators.
	The tuple ``tt`` is formed by the variables of a Laurent polynomial ring.
	
	INPUT:
	
	- ``matriz`` -- An $n\\times n$ invertible matrix over $\\mathbb{Z}$.
	- ``smith`` -- An $m\\times n$ diagonal matrix over $\\mathbb{Z}$.
	- ``tt`` -- Tuple of variables of a Laurent polynomial ring.
	
	
	OUTPUT:

	A dictionnary realizing the relations of the variables and an ideal generated by $t_i^{d_i}-1$ for the
	diagonal terms $>1$.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: R=LaurentPolynomialRing(QQ,'t',4)
		sage: A=random_matrix(ZZ,3,4)
		sage: S,U,V=A.smith_form()
		sage: S
		[1 0 0 0]
		[0 1 0 0]
		[0 0 2 0]
		sage: CambioVarSmith(V,S,R.gens())
		({t3: t3^357, t2: t2*t3^-86, t1: t3^53, t0: t2*t3^-2043}, [t2^2 - 1])

	"""
	A=matriz
	Sm=smith
	n=A.nrows()
	n1=min(Sm.dimensions())
	I=[]
	for i in range(n):
		if i<n1 and Sm[i,i]>1:
			I.append(tt[i]^Sm[i,i]-1)
	dic={}
	for i in range(n):
		T=1
		for j in range(n):
			if j<n1 and Sm[j,j]>0:
				m=A[i,j].quo_rem(Sm[j,j])[1]
			else:
				m=A[i,j]    
			T=T*(tt[j]^m)
		dic.update({tt[i]:T})
	return (dic,I)

def caracter(x,R,cambio):
	r"""
	Given an element `x` in the group algebra of a free group, we apply to it a character
	defined by a dictionnary ``cambio`` to obtain  an element in a Laurent polynomial ring `R`.
	
	INPUT:
	
	- ``x`` -- An element in $\\mathbb{Z}[\\mathbb{F}_n]$.
	- ``R`` -- A ring, usually polynomial or Laurent polynomial.
	- ``cambio`` -- A dictionnary changing the variables.
	
	
	OUTPUT:

	The image by the character of the element $x$.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(2)
		sage: A=F.algebra(QQ)
		sage: a=A(F([1,2,-1,-2]))+3*A(F([2,1,2]))
		sage: R.<t0,t1>=LaurentPolynomialRing(QQ)
		sage: cambio={t1:t0,t0:t0^-2*t1}
		sage: caracter(a,R,cambio)
		3*t1 + 1

	"""
	tt=R.gens()
	res=R(0)
	for i in list(x):
		res+=i[1]*prod([tt[abs(j)-1]^sign(j) for j in i[0].Tietze()])
	return (res.subs(cambio))

def unidades(f,R):
	r"""
	If a Laurent polynomial  ``f`` in a ring ``R`` is non zero, it can be written as a unit and 
	a polynomial divided by no variable.
	
	INPUT:
	
	- ``f`` -- A Laurent polynomial.
	- ``R`` -- The ring of $f$.
	
	
	OUTPUT:

	Either zero or the maximal unit.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: R=LaurentPolynomialRing(QQ,3,'t')
		sage: f=random_vector(QQ,4)
		sage: g1=random_vector(ZZ,4)
		sage: g2
		(0, 0, 2, 1)
		sage: R=LaurentPolynomialRing(QQ,3,'t')
		sage: f=random_vector(QQ,4)
		sage: g1=random_vector(ZZ,4)
		sage: g2=random_vector(ZZ,4)
		sage: p=sum(i*R.gen(0)^j*R.gen(1)^k for (i,j,k) in zip(f,g1,g2))
		sage: p
		-3*t0^3*t1^8 - t0*t1^-1 - 9/2*t1^-13
		sage: unidades(p,R)
		t1^-13

	"""
	if f==f.parent(0):
		return (0)
	ex=f.exponents()
	n=f.parent().ngens()
	u=1
	for j in range(n):
		mx=min([v[j] for v in ex])
		u=u*R.gen(j)^mx
	return (u)

def unidadeslista(lista,R):
	r"""
	Given a list ``lista`` of Laurent polynomials in a ring ``R``, it extracts the maximal commun unit. To be used with 
	base ring of characteristic zero.
	
	INPUT:
	
	- ``lista`` -- A list of Laurent polynomials.
	- ``R`` -- The ring of the elements in ``lista``.
	
	
	OUTPUT:

	Either zero or the maximal unit.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: R=LaurentPolynomialRing(QQ,3,'t')
		sage: F=[random_vector(QQ,4) for i in range(4)]
		sage: G1=[random_vector(ZZ,4) for i in range(4)]
		sage: G2=[random_vector(ZZ,4) for i in range(4)]
		sage: P=[sum(i*R.gen(0)^j*R.gen(1)^k for (i,j,k) in zip(f,g1,g2)) for (f,g1,g2) in zip(F,G1,G2)]
		sage: P
		[t0^5*t1^-1 - t0^2 + 9*t0*t1 + 1/5*t0^-6*t1^-1,
		3/8*t0^-1*t1^2 + 3*t0^-1*t1^-1 - 1/3*t0^-3*t1,
		-2*t0*t1^3 + 16*t1^4,
		-207/704*t1^-1 - 2*t0^-2]
		sage: unidadeslista(P,R)
		t0^-6*t1^-1

	"""
	monomios=sum([unidades(_,R) for _ in lista])
	return (unidades(monomios,R))

def unidadesmatriz(matriz,R,O='R'):
	r"""
	This function extracts the maximal commun unit for each row (or column) of a matrix ``matriz`` 
	of polynomials or Laurent polynomials in a ring ``R``. The default choice of ``O`` is ``R`` (rows); the column 
	choice follow if ``O='C'``.
	
	INPUT:
	
	- ``matrix`` -- A matrix of Laurent polynomials.
	- ``R`` -- The ring of the elements in ``lista``.
	- ``O`` -- A chain of one character, ``'R'`` (for rows), ``'C'`` (for columns) (default: 'R').
	
	OUTPUT:

	A matrix of polynomials with no variable as common divisor in the rows or columns.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: R=LaurentPolynomialRing(QQ,3,'t')
		sage: F=[random_vector(QQ,4) for i in range(6)]
		sage: G1=[random_vector(ZZ,4) for i in range(6)]
		sage: G2=[random_vector(ZZ,4) for i in range(6)]
		sage: P=[sum(i*R.gen(0)^j*R.gen(1)^k for (i,j,k) in zip(f,g1,g2)) for (f,g1,g2) in zip(F,G1,G2)]
		sage: A=Matrix(R,3,P)
		sage: A
		[       39/2*t0^3*t1^-1 - 1/23*t0*t1^-1 - 12*t0^-2*t1 -3/2*t1^46 + 48*t0^4*t1^2 + 2*t0^-1*t1 + t0^-1*t1^-3]
		[              -t0^10*t1^5 - 4*t0^6 - 4/221*t0*t1^-20                          2*t0*t1^-1 - 6 + 1/63*t0^-2]
		[                                            49*t1^-1     -1/2 - 1/3*t0^-3*t1^3 - 3/8*t1^-1 + 9*t0^2*t1^-5]
		sage: unidadesmatriz(unidadesmatriz(A,R),R,O='C')
		[         39/2*t0^5*t1^2 - 1/23*t0^3*t1^2 - 12*t1^4    -3/2*t0^2*t1^49 + 48*t0^6*t1^5 + 2*t0*t1^4 + t0]
		[          -t0^12*t1^25 - 4*t0^8*t1^20 - 4/221*t0^3           2*t0^3*t1^19 - 6*t0^2*t1^20 + 1/63*t1^20]
		[                                      49*t0^3*t1^4 -1/2*t0^3*t1^5 - 1/3*t1^8 - 3/8*t0^3*t1^4 + 9*t0^5]

	"""
	if O=='C':
		return (unidadesmatriz(matriz.transpose(),R).transpose())
	A=matriz.change_ring(R)
	n=A.nrows()
	L=A.rows()
	U=[unidadeslista(_,R) for _ in L]
	for i in range(n):
		if U[i]!=A.base_ring()(0):
			L[i]=U[i]^-1*L[i]
	return (Matrix(L))




def Hay_unidades(A,R,S,dividir=False):
	r"""
	This function checks if there is a unit in a matrix of polynomials or Laurent polynomials in a ring ``R``. There is an optional parameter ``dividir`` having as default value ``False``.
	
	INPUT:
	
	- ``A`` -- A matrix of polynomials or Laurent polynomials.
	- ``R`` -- A ring of Laurent polynomials.
	- ``S`` -- The associated ring of polynomials.
	- ``dividir`` -- A boolean.
	
	OUTPUT:

	If there is no unit (over the integers if ``dividir`` is ``False``) the output is ``None``; if not, the position of the first unit.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: R=LaurentPolynomialRing(QQ,1,'t')
		sage: A=Matrix(R,2,2,[R.gen(0)^-1,2,R.gen(0),R.gen(0)-1])
		sage: Hay_unidades(A,R,R.polynomial_ring())
			(0,0)                
		sage: A.rescale_col(0,2*R.gen(0)^2)
		sage: Hay_unidades(A,R,R.polynomial_ring())==None
			(0,0)                

	"""
	res=None
	n,m=A.dimensions()
	i=0
	j=0
	while res==None and j<m:
		if R(A[i,j]).is_unit():
			if not dividir and A[i,j].coefficients()[0]^2==1:
				res=(i,j)
			elif dividir:
				res=(i,j)
		i+=1
		if i==n:
			j+=1
			i=0
	return (res)


def reducir_matriz(A,R,S,ideal,dividir=False):
	r"""This function performs standard matrix operations eliminating such that the new matrix represents the same module

	INPUT:

	- ``matrix`` -- A matrix of Laurent polynomials.
	- ``R`` -- Laurent polynomial ring.
	- ``S`` -- Polynomial ring associated to ``R``.
	- ``ideal`` Torsion ideal in ``S``.
	- ``dividir`` -- A boolean.

	OUTPUT:

	The reduced matrix.

	EXAMPLES:

		sage: R=LaurentPolynomialRing(QQ,1,'t')
		sage: A=Matrix(R,2,2,[R.gen(0)^-1,2,R.gen(0),R.gen(0)-1])
		sage: R=LaurentPolynomialRing(QQ,1,'t')
		sage: S=R.polynomial_ring()
		sage: reducir_matriz(A,R,S,S.ideal([0]))
		[-2*t^2 + t - 1]

	"""
	n,m=A.dimensions()
	A1=unidadesmatriz(unidadesmatriz(A,R),R,O='C').change_ring(S)
	A1=Matrix(A1.nrows(),[v.reduce(S.ideal(ideal)) for v in A1.list()])
	A1=A1.change_ring(R)
	res=Hay_unidades(A1,R,S,dividir)
	while res!=None:
		i,j=res
		A1.swap_rows(i,0)
		A1.swap_columns(j,0)
		a=[i for i in range(1,m) if A1[0,i]==0]
		b=[i for i in range(1,n) if A1[i,0]==0]
		if a<=b:
			for j in range(1,m):
				A1.add_multiple_of_column(j,0,-A1[0,j]*A1[0,0]^-1)
		else:
			for i in range(1,n):
				A1.add_multiple_of_row(i,0,-A1[i,0]*A1[0,0]^-1)
		A1=A1.delete_rows([0]).delete_columns([0])
		A1=unidadesmatriz(unidadesmatriz(A1,R),R,O='C').change_ring(S)
		A1=Matrix(A1.nrows(),[v.reduce(S.ideal(ideal)) for v in A1.list()])
		A1=A1.change_ring(R)
		res=Hay_unidades(A1,R,S,dividir)
		n,m=A1.dimensions()
	return (A1)

def grafoplumbing(grafo,selfint,generos=None,flechas=None):
	r"""
	Given a (simplicial) graph with two weights (self-intersection and genus) and two types of vertices (standard ones and arrows), this function provides a finite presentation of the fundamental group of the complement of the link (eventually void) of a graph manifold represented by this plumbin graph.
	
	INPUT:
	
	- ``grafo`` -- A simplicial graph.
	- ``selfint`` -- A list of length the number $n$ of vertices of the graph, containing the self-intersections of those vertices.
	- ``generos`` -- A list of $n$ non-negative integers for the genera: the value ``None`` means all of them are zero.
	- ``flechas`` -- A list of $n$ non-negative integers for the number of arrows attached at it vertex: the value ``None`` means all of them are zero.
	
	OUTPUT:

	A presentation of the group


	EXAMPLES:

	This example illustrates a simple use of this function

	::

		sage: gr=graphs.CycleGraph(3)
		sage: inter=[-2,-2,-3]
		sage: g=grafoplumbing(gr,inter)
		Definida negativa:  True
		sage: g
		Finitely presented group < a0, a1, a2, b | a0^-2*b*a1*b^-1*a2, b^-1*a0*b*a1^-2*a2, a0*a1*a2^-3, a0*b*a1*b^-1*a0^-1*b*a1^-1*b^-1 >

	"""

	n=grafo.num_verts()
	if generos==None:
		generos=n*[0]
	if flechas==None:
		flechas=n*[0]
	A=grafo.adjacency_matrix()
	for i in range(n):
		A[i,i]=selfint[i]
	print("Definida negativa: ",(-A).is_positive_definite())
	T=grafo.spanning_trees()[0]
	L=[_ for _ in grafo.edges(labels=False) if _ not in T.edges(labels=False)]
	m=len(L)
	s=sum(generos)
	r=sum(flechas)
	F0=FreeGroup(n,'a')
	F1=FreeGroup(m,'b')
	F2=FreeGroup(s,'c')
	F3=FreeGroup(s,'d')
	F4=FreeGroup(r,'e')
	F=FreeGroup(F0.gens()+F1.gens()+F2.gens()+F3.gens()+F4.gens())
	gensgeneros1=[]
	gensgeneros2=[]
	gensflechas=[]
	inicio=n+m
	for j in range(n):
		final=inicio+generos[j]
		gn1=[F.gen(j) for j in range(inicio,final)]
		gn2=[F.gen(j) for j in range(inicio+s,final+s)]
		inicio=final
		gensgeneros1.append(gn1)
		gensgeneros2.append(gn2)
	inicio=n+m+2*s
	for j in range(n):
		final=inicio+flechas[j]
		gn=[F.gen(j) for j in range(inicio,final)]
		inicio=final
		gensflechas.append(gn)
	F.inject_variables(verbose=False)
	a=F.gens()
	rlk=[]
	for i in range(n):
		w=F(1)
		for j in range(n):
			if (i,j) in L:
				k=L.index((i,j))
				w=w*a[n+k]*a[j]^A[i,j]/a[n+k]
			elif (j,i) in L:
				k=L.index((j,i))
				w=w/a[n+k]*a[j]^A[i,j]*a[n+k]
			else:
				w=w*a[j]^A[i,j]
		for (x,y) in zip(gensgeneros1[i],gensgeneros2[i]):
			w=w*x*y/x/y
		for x in gensflechas[i]:
			w=w*x
		rlk.append(w)
	for (i,j) in L:
		k=L.index((i,j))
		w=a[n+k]*a[j]/a[n+k]
		rlk.append(a[i]*w/a[i]/w)
	for j in range(n):
		for x in gensgeneros1[j]+gensgeneros2[j]:
			rlk.append(a[j]*x/a[j]/x)
	for j in range(n):
		for x in gensflechas[j]:
			rlk.append(a[j]*x/a[j]/x)
	g=F/rlk
	return (g)
	
def AnilloMatriz(grupo):
	r"""
	The Alexander matrix of a group ``grupo`` related to the derived subgroup.
	
	INPUT:
	
	- ``grupo`` -- A finitely presented group.
	
	OUTPUT:

	It prints the dictionnary with the translations of the group generators to the Laurent polynomial variables
	and returns the Laurent polynomial ring, the associate polynomial ring, the 
	Alexander matrix of a group in the group algebra of the abelianization and an ideal which is non-zero
	if the homology has non trivial torsion


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(4)
		sage: AnilloMatriz(B)
		{t2: t2, t0: t2, t1: t2}
		(
		Multivariate Laurent Polynomial Ring in t0, t1, t2 over Rational Field,

		Multivariate Polynomial Ring in t0, t1, t2 over Rational Field,

		[ t2^2 - t2 + 1 -t2^2 + t2 - 1              0]    
		[       -t2 + 1              0         t2 - 1]    
		[             0  t2^2 - t2 + 1 -t2^2 + t2 - 1], []
		)
		sage: F=FreeGroup(2); rlk=[F.gen(0)^2,F.gen(1)^3]; G=F/rlk
		sage: AnilloMatriz(G)
		{t0: t1^3, t1: t1^2}
		(
		Multivariate Laurent Polynomial Ring in t0, t1 over Rational Field,

		Multivariate Polynomial Ring in t0, t1 over Rational Field,

		[       t1^3 + 1               0]            
		[              0 t1^4 + t1^2 + 1], [t1^6 - 1]
		)

	"""
	MR=MatrizAbel(grupo)
	Sm,U,V=MR.smith_form()
	M=grupo.alexander_matrix()
	R=LaurentPolynomialRing(QQ,len(grupo.gens()),'t')
	S=R.polynomial_ring()
	R.inject_variables(verbose=False)
	tt=R.gens()
	cambio,idl=CambioVarSmith(V,Sm,tt)
	print (cambio)
	A=matrix(M.nrows(),M.ncols(),[caracter(_,R,cambio) for _ in M.list()])
	A1=unidadesmatriz(unidadesmatriz(A,R),R,O='C').change_ring(S)
	A1=Matrix(A1.nrows(),[v.reduce(S.ideal(idl)) for v in A1.list()])
	vars=[]
	for a in cambio.keys():
		vr=list(cambio[a].variables())
		vars+=vr
	vars=tuple(Set(vars))
	R1=LaurentPolynomialRing(QQ,vars)
	S1=R1.polynomial_ring()
	idl1=[S1(S(_)) for _ in idl]
	A2=A1.change_ring(S1)
	return (R1,S1,A2,idl1)


def pseudo_coxeter(lista):
	r"""
	Given a list ``lista`` of braids, it returns its pseudo Coxeter element. In this convention, it is the reversed
	product.
	
	INPUT:
	
	- ``lista`` -- A list of braids.
	
	OUTPUT:

	The reversed product.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(5)
		sage: pseudo_coxeter(B.gens())
		s3*s2*s1*s0

	"""
	return (prod([_ for _ in reversed(lista)]))

def gapperm(listatupla):
	r"""
	Given a list ``listatupla`` of tuples of positive integers (or just one tuple), it returns the corresponding
	GAP permutation. The integers may be repeted
	
	INPUT:
	
	- ``listatupla`` -- A list of tuples of positive integers, or a tuple.
	
	OUTPUT:

	A GAP permutation.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: gapperm([(2,1),(1,5)])
		(1,2,5)
		sage: gapperm((2,1))
		(1,2)

	"""
	if type(listatupla)==tuple:
		return (gap(PermutationGroupElement(listatupla)))
	res=gap(PermutationGroupElement(()))
	for tupla in listatupla:
		res=res*gap(PermutationGroupElement(tupla))
	return (res)

def generadores_trenzas(lista,hilos):
	r"""
	Given a list ``lista`` of list or tuples encoding permutations in $\\Sigma_n$, $n=$ ``hilos``, it returns
	a list of generators of the preimage of the subgroup generated by these generators under the natural
	map $\\mathbb{F}_n\to\\Sigma_n$
	
	INPUT:
	
	- ``lista`` -- A list of of list  of tuples of positive integers (or tuples).
	- ``hilos`` -- Number of strands of the braid group
	
	OUTPUT:

	Generators of the subgroup of the braid group realizing the permutations in the list.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: generadores_trenzas([()],2)
		[s^-2]
		sage: generadores_trenzas([(1,2)],3)
		[s0, s1^-2]

	"""
	B0=BraidGroup(hilos)
	Bg=B0.gap()
	Sg=gap(hilos).SymmetricGroup()
	genB=Bg.GeneratorsOfGroup()
	im=[]
	for j in range(hilos-1):
		im.append(gapperm((j+1,j+2)))
	homg=Bg.GroupHomomorphismByImages(Sg,genB,im)
	k0=Sg.Subgroup(map(gapperm,lista))
	genS=homg.PreImage(k0).GeneratorsOfGroup()
	generadores=[B0(_.UnderlyingElement().TietzeWordAbstractWord().sage()) for _ in genS]
	return (generadores)

def GapConvert(objeto):
	r"""
	Convert ``object`` to GAP.
	
	INPUT:
	
	- ``objeto`` -- Some group related object.
	
	OUTPUT:

	The object ``objeto`` in GAP.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: GapConvert(BraidGroup(3))
		<fp group of size infinity on the generators [ s0, s1 ]>
		
	"""
	tipo=str(type(objeto)).upper()
	if 'GAP' in tipo and 'LINEAR' not in tipo:
		return (objeto)
	elif 'PERM' in tipo:
		return (gap(objeto))
	else:
		return (objeto.gap())
	


def HomomorphismGroups(G1,G2,im):
	r"""
	Given two (SAGEMATH or GAP) groups ``G1``, ``G2`` (with generator systems) and a lista ``im`` of elements
	in ``G2`` (as many as the generators) of ``G1``, it constructs the associated GAP homomorphism.
	
	INPUT:
	
	- ``G1`` -- SAGEMATH or GAP finitely presented group.
	- ``G2`` -- SAGEMATH or GAP finitely presented group.
	- ``im`` -- list of elements in ``G2``
	
	OUTPUT:

	GAP homomorphism.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: HomomorphismGroups(B.gap(),gap(g),[gapperm(_) for _ in [(1,2),(2,3)]])
		[ s0, s1 ] -> [ (1,2), (2,3) ]
		sage: F=FreeGroup(2);B=BraidGroup(3)
		sage: HomomorphismGroups(F,B,B.gens())
		[ x0, x1 ] -> [ s0, s1 ]
		sage: HomomorphismGroups(F.gap(),B.gap(),B.gap().GeneratorsOfGroup())
		[ x0, x1 ] -> [ s0, s1 ]
		sage: HomomorphismGroups(F,B.gap(),B.gap().GeneratorsOfGroup())
		[ x0, x1 ] -> [ s0, s1 ]
		sage: g=SymmetricGroup(3)
		sage: HomomorphismGroups(B,g,[PermutationGroupElement(_) for _ in [(1,2),(2,3)]])
		[ s0, s1 ] -> [ (1,2), (2,3) ]
		sage: HomomorphismGroups(B.gap(),gap(g),[gapperm(_) for _ in [(1,2),(2,3)]])
		[ s0, s1 ] -> [ (1,2), (2,3) ]
		
	"""
	g1=G1.gap()
	g2=G2.gap()
	gen=g1.GeneratorsOfGroup()
	im0=[_.gap() for _ in im]
	hom=g1.GroupHomomorphismByImages(g2,gen,im0)
	return (hom)


def HomImageEl(hom,el):
	r"""
	It returns the image in the target of an element ``el`` in the source of a GAP group homomorphism ``hom``.
	
	INPUT:
	
	- ``hom`` -- GAP group homomorphism.
	- ``el`` -- element in the source of ``hom``
	
	OUTPUT:

	The image of ``el`` by ``hom``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(2);B=BraidGroup(3)
		sage: hom=HomomorphismGroups(F,B,B.gens())
		sage: HomImageEl(hom,F([1,-2]))
		s0*s1^-1
		
	"""
	return (hom.Image(GapConvert(el)))

def HomImageSub(hom,lista):
	r"""
	Given a GAP group homomorphism ``hom`` and a list ``lista`` of elements in the source, we provide the GAP subgroup image of
	the subgroup generated by ``lista``.
	
	INPUT:
	
	- ``hom`` -- GAP group homomorphism.
	- ``lista`` -- elements in the source of ``hom``
	
	OUTPUT:

	The image of the subgroup generated by ``lista`` by ``hom``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: F=FreeGroup(3);B=BraidGroup(4)
		sage: hom=HomomorphismGroups(F,B,B.gens())
		sage: HomImageSub(hom,[F.gen(0)^2,F.gen(1)/F.gen(2)])
		Group([ s0^2, s1*s2^-1 ])
		
	"""
	listag=map(GapConvert,lista)
	listaim=[hom.Image(_) for _ in listag]
	gs=hom.Range()
	im=gs.Subgroup(listaim)
	return (im)


def Centralizador(ggap,elgap):
	r"""
	Given a GAP group ``ggap`` and an element ``elgap`` (maybe not in ``ggap`` but in some parent group),
	this function returns the centralizer of the element in the group.
	
	INPUT:
	
	- ``ggap`` -- GAP group .
	- ``elgap`` -- element in some parent group
	
	OUTPUT:

	Centralizer of ``elgap`` in ``ggap``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: g=gap(DihedralGroup(4))
		sage: a=gapperm((1,2))
		sage: a in g
		False
		sage: Centralizador(g,a)
		Group( [ (1,2)(3,4) ] )
		
	"""
	return (ggap.Centralizer(elgap))

def Normalizador(ggap,subgap):
	r"""
	Given a GAP group ``ggap`` and an group ``subgap`` (maybe not in ``ggap`` but in some parent group),
	this function returns the normalizer of ``subgap`` in ``ggap``.
	
	INPUT:
	
	- ``ggap`` -- GAP group.
	- ``subgap`` -- GAP group with some common parent with ``ggap``
	
	OUTPUT:

	Normalizer of ``subgap`` in ``ggap``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: g=gap(DihedralGroup(4))
		sage: s=gap(SymmetricGroup(4))
		sage: Normalizador(s,g)
		Group( [ (2,4), (1,4,3,2), (1,3)(2,4) ] )
		
	"""
	return (ggap.Normalizer(subgap))

def Hurwitz0(B,gensS,br1,br2,GS,lista):
	r"""
	The group ``B`` is a braid group, ``gensS`` is a system of generators of the braid subgroup which fixes some
	subsets of strands. The lists ``br1`` and ``br2`` are braid monodromies (lists of braids in ``B``,
	``GS`` is a finite group and ``lista`` is a family of elements of ``GS`` defining a representation
	of ``B`` into ``GS``.
	
	INPUT:
	
	- ``B`` -- A braid group.
	- ``gensS`` -- a list of elements in ``B``
	- ``br1`` -- a list of $r$ braids in ``B``
	- ``br2`` -- another list of $r$ braids in ``B``
	- ``GS`` -- a finite GAP or SAGEMATH group 
	- ``lista`` -- a list of elements of ``GS``
	
	OUTPUT:

	If the images of the braid monodromies by the morphism cannot have either the
	same pseudo Coxeter element $c$ or the same monodromy group $G$, up to conjugation, it returns ``None``.
	If it is the case, a list with the following elements: 
	the image of ``br1`` by the morphism, a list with all the conjugated of the image of ``br2`` 
	with the same  pseudo Coxeter element and the same monodromy group as the one of ``br1``,
	and the group that preserves $c$ and $G$. Finally it returns also the image representation group as a GAP permutation group
	and the list of permutations defining the representation.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(3)
		sage: gensS=[B.gen(0)]
		sage: br1=[B([1,2])^2,B([1])]; br2=[B([1,2])^2,B([2])]
		sage: GS=SymmetricGroup(3); 
		sage: lista=[PermutationGroupElement(_) for _ in [(1,2),(2,3)]]
		sage: Hurwitz0(B,gensS,br1,br2,GS,lista)
		Coxeter distinto
		sage: br2=[B([2,1])^2,B([2])]
		sage: Hurwitz0(B,gensS,br1,br2,GS,lista)
		[[(1,2,3), (1,2)],
		[[(1,3,2), (2,3)]],
		Group(()),
		Group( [ (1,2), (2,3) ] ),
		[(1,2), (2,3)]]
		
	"""
	cx1,cx2=map(pseudo_coxeter,[br1,br2])
	hom=HomomorphismGroups(B,GS,lista)
	listagap=map(GapConvert,lista)
	grupo=GapConvert(GS).Subgroup(listagap)
	GSg=HomImageSub(hom,gensS)
	cxg1,cxg2=[hom.Image(GapConvert(_)) for _ in [cx1,cx2]]
	brg1,brg2=[[hom.Image(GapConvert(_)) for _ in br] for br in [br1,br2]]
	M1g=HomImageSub(hom,br1)
	coxeter=GSg.IsConjugate(cxg2,cxg1).sage()
	#CSg=GSg.Centralizer(cxg1)
	CSg=Centralizador(GSg,cxg1)
	if coxeter:
		vale=GSg.RepresentativeAction(cxg2,cxg1)
		segundo0=[vale^-1*_*vale for _ in brg2]
		M2g=hom.Range().Subgroup(segundo0)
		normales=CSg.IsConjugate(M1g,M2g)
		if not bool(normales):
			return (None)
		vale=CSg.RepresentativeAction(M2g,M1g)
		#NSg=CSg.Normalizer(M1g)
		NSg=Normalizador(CSg,M1g)
		segundo1=[vale^-1*_*vale for _ in segundo0]
		#M2g=hom.Range().Subgroup(segundo1)
		segundo=[[g0*_/g0 for _ in segundo1] for g0 in NSg.Elements()]
		#return ([bool(prod([cxg1==prod([h for h in reversed(_)]) for _ in segundo])),brg1,segundo,NSg])
		return ([brg1,segundo,NSg,grupo,listagap])
	print ('Coxeter distinto')
	return (None)

def HurwitzTrenza(trenza,lista,grupo):
	r"""
	A braid ``trenza`` with $r$ strands acts by Hurwitz moves on a list ``lista`` of $r$ elements of a group ``G``.
	
	INPUT:
	
	- ``trenza`` -- A braid with $r$ strands.
	- ``lista`` -- a list of $r$ elements of the group ``G``
	- ``G`` -- a GAP group

	
	OUTPUT:

	If the images of the braid monodromies by the morphism cannot have either the
	same pseudo Coxeter element $c$ or the same monodromy group $G$, up to conjugation, it returns ``None``.
	If it is the case,a list with the following elements: 
	the image of ``br1`` by the morphism, a list with all the conjugated of the image of ``br2`` 
	with the same  pseudo Coxeter element and the same monodromy group as the one of ``br1``,
	and the group that preserves $c$ and $G$.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B3=BraidGroup(3)
		sage: B4=BraidGroup(4)
		sage: L=list(map(GapConvert,B4.gens()))
		sage: HurwitzTrenza(B3.gen(0),L,B4.gap())
		(s1, s1*s0*s1^-1, s2)
		
	"""
	k=trenza.parent().strands()
	fgens=FreeGroup(k).gens()
	fims=[LibreTrenza(_,trenza).Tietze() for _ in fgens]
	res=[]
	for tupla in fims:
		res0=grupo.One()
		for j in tupla:
			res0=res0*lista[j.abs()-1]^j.sign()
		res.append(res0)
	return (tuple(res))


def Hurwitz1(B,gensS,generadores,br1,br2,GS,lista):
	r"""
	The group ``B`` is a braid group, ``gensS`` is a system of generators of the braid subgroup which fixes some
	subsets of strands. The list ``generadores`` is a generator of the subgroup of allowed Hurwitz moves
	while the lists ``br1`` and ``br2`` are braid monodromies (lists of braids in ``B``,
	``GS`` is a finite group and ``lista`` is a family of elements of ``GS`` defining a representation
	of ``B`` into ``GS``.
	
	INPUT:
	
	- ``B`` -- A braid group.
	- ``gensS`` -- a list of elements in ``B``
	- ``generadores`` -- a list of elements in the braid group of $r$ strands
	- ``br1`` -- a list of $r$ braids in ``B``
	- ``br2`` -- another list of $r$ braids in ``B``
	- ``GS`` -- a a finite GAP or SAGEMATH group
	- ``lista`` -- a list of elements of ``GS``
	
	OUTPUT:

	It constructs the Hurwitz orbit of the images of ``br1`` in $G:=$ ``GS``. If the image of ``br2`` appears in the orbit,
	then it returns ``None``; if not it returns the orbit in $G^r/G$ (quotient by conjugation).
	Finally it returns also the image representation group as a GAP permutation group
	and the list of permutations defining the representation.



	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(3)
		sage: gensS=[B.gen(0)]
		sage: br1=[B([1,2])^2,B([1]),B([2])]
		sage: br2=[B([2,1])^2,B([2]),B([1])]
		sage: GS=SymmetricGroup(3)
		sage: lista=[PermutationGroupElement(_) for _ in [(1,2),(2,3)]]
		sage: generadores=generadores_trenzas([(1,2)],3)
		sage: Hurwitz1(B,gensS,generadores,br1,br2,GS,lista)
		aumentar  1
		nuevos  2
		orbita  3
		aumentar  2
		Orbitas coinciden
		[((1,2,3), (1,3), (1,2))]
		
	"""
	hurwitz=Hurwitz0(B,gensS,br1,br2,GS,lista)
	if hurwitz==None:
		print ("Orbitas no coinciden")
		return (None)
	brg1,segundo,gr,grupo,listagap=hurwitz
	segundo=map(tuple,segundo)
	orbita=[tuple(brg1)]
	aumentar=[tuple(brg1)]
	if brg1 in segundo:
		print ("Orbitas coinciden")
		return (segundo)
	control0=True
	while control0:
		nuevos=[]
		print ("aumentar ",len(aumentar))
		for elt in aumentar:
			for kj in generadores:
				elt1=HurwitzTrenza(kj,elt,GapConvert(GS))
				control=elt1 in segundo
				if control:
					print ("Orbitas coinciden")
					return (segundo)
				orbcnj=[tuple([g^-1*_*g for _ in elt1]) for g in gr.Elements()]
				#control1=len(Set(orbita).intersection(Set(orbcnj)))==0
				control1=True
				c=0
				while control1 and c<len(orbcnj):
					_=orbcnj[c]
					control1=control1 and _ not in orbita
					c+=1
				if control1:
					orbita.append(elt1)
					nuevos.append(elt1)    
			#print ("nuevos ",len(nuevos),"\n lugar ",aumentar.index(elt),"\n")
		print ("nuevos ",len(nuevos))
		print ("orbita ",len(orbita))
		aumentar=copy(nuevos)
		control0=len(nuevos)>0  
	print ("Orbitas no coinciden")
	return (orbita,segundo,gr,grupo,listagap)


def HurwitzOrbita(B,gensS,generadores,br1,GS,lista):
	r"""
	The group ``B`` is a braid group, ``gensS`` is a system of generators of the braid subgroup which fixes some
	subsets of strands. The list ``generadores`` is a generator of the subgroup of allowed Hurwitz moves
	while the list ``br1`` is a braid monodromy (lists of braids in ``B``,
	``GS`` is a finite group and ``lista`` is a family of elements of ``GS`` defining a representation
	of ``B`` into ``GS``.
	
	INPUT:
	
	- ``B`` -- A braid group.
	- ``gensS`` -- a list of elements in ``B``
	- ``generadores`` -- a list of elements in the braid group of $r$ strands
	- ``br1`` -- a list of $r$ braids in ``B``
	- ``GS`` -- a finite GAP or SAGEMATH group
	- ``lista`` -- a list of elements of ``GS``
	
	OUTPUT:

	It constructs the Hurwitz orbit of the images of ``br1`` in $G:=$ ``GS``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(3)
		sage: gensS=[B.gen(0)]
		sage: br1=[B([1,2])^2,B([1]),B([2])]
		sage: GS=SymmetricGroup(3)
		sage: lista=[PermutationGroupElement(_) for _ in [(1,2),(2,3)]]
		sage: generadores=generadores_trenzas([(1,2)],3)
		sage: HurwitzOrbita(B,gensS,generadores,br1,GS,lista)
		aumentar  1
		nuevos  2
		orbita  3
		aumentar  2
		nuevos  4
		orbita  7
		aumentar  4
		nuevos  3
		orbita  10
		aumentar  3
		nuevos  2
		orbita  12
		aumentar  2
		nuevos  0
		orbita  12
		([((1,2,3), (1,2), (2,3)),
		((1,2), (1,3,2), (2,3)),
		((1,2,3), (2,3), (1,3)),
		((1,3,2), (2,3), (2,3)),
		((1,2), (1,2,3), (1,2)),
		((2,3), (1,3,2), (1,3)),
		((1,2,3), (1,3), (1,2)),
		((2,3), (1,2,3), (2,3)),
		((1,3,2), (1,3), (1,3)),
		((1,3), (1,3,2), (1,2)),
		((1,3), (1,2,3), (1,3)),
		((1,3,2), (1,2), (1,2))],
		Group(()),
		Group( [ (1,2), (2,3) ] ),
		[(1,2), (2,3)])
		
	"""
	cx1=pseudo_coxeter(br1)
	hom=HomomorphismGroups(B,GS,lista)
	listagap=map(GapConvert,lista)
	grupo=GapConvert(GS).Subgroup(listagap)
	GSg=HomImageSub(hom,gensS)
	cxg1=hom.Image(GapConvert(cx1))
	brg1=[hom.Image(GapConvert(_)) for _ in br1]
	M1g=HomImageSub(hom,br1)
	CSg=Centralizador(GSg,cxg1)
	gr=Normalizador(CSg,M1g)
	orbita=[tuple(brg1)]
	aumentar=[tuple(brg1)]
	control0=True
	while control0:
		nuevos=[]
		print ("aumentar ",len(aumentar))
		iterador=xmrange_iter([aumentar,generadores])
		for (elt,kj) in iterador:
			elt1=HurwitzTrenza(kj,elt,GapConvert(GS))
			orbcnj=[tuple([g^-1*_*g for _ in elt1]) for g in gr.Elements()]
			control1=True
			c=0
			while control1 and c<len(orbcnj):
				_=orbcnj[c]
				control1=control1 and _ not in orbita
				c+=1
			if control1:
				orbita.append(elt1)
				nuevos.append(elt1)    
		print ("nuevos ",len(nuevos))
		print ("orbita ",len(orbita))
		aumentar=copy(nuevos)
		control0=len(nuevos)>0  
	return (orbita,gr,grupo,listagap)
	
	
def burauperm(nm,m,hilos,field=False):
	r"""
	From Burau representation it is possible to construct a finite representation
	of the braid group `B` in ``hilos`` strands. It is done using $\\mathbb{Z}/n_m$, $n_m=$ ``nm`` and $t\\equiv m$, if ``field`` is False,
	or $\\mathbb{F}_{n_m}$ and $t=\\varphi^m$ ($\\varphi$ a generator of the field for $m$ a non trivial 
	prime power), if field is True.
	
	INPUT:
	
	- ``nm`` -- A natural number or a non prime prime power.
	- ``m`` -- an integer
	- ``hilos`` -- number of strands
	- ``field`` -- a boolean (default: False)
	
	OUTPUT:

	A permutation group and a list of $h-1$ elements of the group, $h=$ ``hilos``.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: burauperm(3,2,3)
		(Group([ (2,5,8)(3,9,6), (2,4,9)(3,7,5) ]), [(2,5,8)(3,9,6), (2,4,9)(3,7,5)])
		sage: burauperm(4,1,3,field=True)
		(Group([ (2,4,3), (1,3,4) ]), [(2,4,3), (1,3,4)])
		
	"""
	if field:
		R0=GF(nm,'p')
		mm=R0.gen()^m
	else:
		R0=Zmod(nm)
		mm=R0(m)
	G=GL(hilos-1,R0)
	L.<t>=LaurentPolynomialRing(ZZ)
	A0=[]
	for s in BraidGroup(hilos).gens():
		A=G(Matrix(hilos-1,[R0(_.subs(t=mm)) for _ in s.burau_matrix(reduced=True).list()]))
		A0.append(A)
	iso1=G.gap().IsomorphismPermGroup()
	g1=iso1.Image()
	lista1=[iso1.Image(GapConvert(_)) for _ in A0]
	g1a=g1.Subgroup(lista1)
	iso2=g1a.SmallerDegreePermutationRepresentation()
	g2=iso2.Image()
	lista2=[iso2.Image(GapConvert(_)) for _ in lista1]
	return (g2,lista2)

def compruebaorbita(orbita,lista):
	r"""
	Given lists ``orbita`` and ``lista`` of tuples of $r$ elements in a group it checks if there are common elements.
	
	INPUT:
	
	- ``orbita`` -- A list of tuples of $r$ elements in a group.
	- ``lista`` -- Another list of tuples of $r$ elements in the same group.
	
	OUTPUT:

	It returns ``False`` if no element in common.


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(3)
		sage: gensS=[B.gen(0)]
		sage: br1=[B([1,2])^2,B([1]),B([2])]
		sage: GS=SymmetricGroup(3)
		sage: lista=[PermutationGroupElement(_) for _ in [(1,2),(2,3)]]
		sage: generadores=generadores_trenzas([(1,2)],3)
		sage: res1=Hurwitz1(B,gensS,generadores,br1,br2,GS,lista)
		aumentar  1
		nuevos  2
		orbita  3
		aumentar  2
		Orbitas coinciden
		sage: res2=HurwitzOrbita(B,gensS,generadores,br1,GS,lista)
		aumentar  1
		nuevos  2
		orbita  3
		aumentar  2
		nuevos  4
		orbita  7
		aumentar  4
		nuevos  3
		orbita  10
		aumentar  3
		nuevos  2
		orbita  12
		aumentar  2
		nuevos  0
		orbita  12
		sage: segundo=res1
		sage: orbita=res2[0]
		sage: compruebaorbita(orbita,segundo)
		True
		
	"""
	res=False
	l=len(lista)
	k=0
	while not res and k<l:
		res=res or lista[k] in orbita
		k=k+1
	return (res)


def orbitapermutacion(orbita,grupo,generadores):
	r"""
	The list ``orbita`` contains tuples of $r$ elements of a group $G$, ``grupo`` is a subgroup $H$ of $G$
	and ``generadores`` is a list of generators of the subgroup of allowed Hurwitz moves.
	
	INPUT:
	
	- ``orbita`` -- tuples of $r$ elements of a group.
	- ``grupo`` -- a subgroup of the above group
	- ``generadores`` -- a list of elements in the braid group of $r$ strands
	
	OUTPUT:

	It returns the permutations defined by ``generadores``by its Hurwitz action on ``orbita``
	as elements in $G^r/H$ (conjugation action).


	EXAMPLES:

	This example illustrates a simple use of this function

	::
	
		sage: B=BraidGroup(3)
		sage: gensS=[B.gen(0)]
		sage: br1=[B([1,2])^2,B([1]),B([2])]
		sage: GS=SymmetricGroup(3)
		sage: lista=[PermutationGroupElement(_) for _ in [(1,2),(2,3)]]
		sage: generadores=generadores_trenzas([(1,2)],3)
		sage: res=HurwitzOrbita(B,gensS,generadores,br1,GS,lista)
		aumentar  1
		nuevos  2
		orbita  3
		aumentar  2
		nuevos  4
		orbita  7
		aumentar  4
		nuevos  3
		orbita  10
		aumentar  3
		nuevos  2
		orbita  12
		aumentar  2
		nuevos  0
		orbita  12
		sage: orbita=res[0]
		sage: grupo=res[2]
		sage: orbitapermutacion(orbita,grupo,generadores)
		[( 1, 2, 4, 5)( 3, 6, 9, 8)( 7,10,12,11), ( 2, 5)( 6, 8)(10,11)]
		
	"""
	grel=grupo.Elements()
	orbitaconj=[tuple(sorted([[g*A/g for A in _] for g in grel])) for _ in orbita]
	permutaciones=[]
	orbitaconjnum=list(enumerate(orbitaconj, 1))
	for s in generadores:
		im=[HurwitzTrenza(s,_,GapConvert(grupo)) for _ in orbita]
		imconj=[tuple(sorted([[g*A/g for A in _] for g in grel])) for _ in im]
		imconjnum=list(enumerate(imconj, 1))
		porbita=Permutation([pair[0] for pair in sorted(orbitaconjnum, key=lambda x: x[1])])
		pim=Permutation([pair[0] for pair in sorted(imconjnum, key=lambda x: x[1])])
		permutaciones.append(gap(pim.inverse()*porbita))
	return (permutaciones)

