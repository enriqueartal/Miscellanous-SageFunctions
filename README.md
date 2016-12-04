# Miscellanous-SageFunctions
Some auxiliar functions in Sagemath
It contains right now a file named MonodromyGroupFunctions.sage with functions mostly related with braid monodromy. Functions are documented. Last update: 11-06-2016
I use the following code to use these functions in a worksheet or a sage session:

  import os
  DIR='path_where_the_file_is_universally_accessible/'
  SAGE_LOCATION='path_where_sagemath_is_installed/'
  if 'MonodromyGroupFunctions.sage.py' in os.listdir(DIR):
      NO=os.system('rm '+DIR+'MonodromyGroupFunctions.sage.py')
  NO=os.system(SAGE_LOCATION+'local/bin/sage-preparse '+DIR+'MonodromyGroupFunctions.sage')
  os.chmod(DIR+'MonodromyGroupFunctions.sage.py',0o666)
  load(DIR+'MonodromyGroupFunctions.sage.py')
 

In a sage worksheet type Name_Of_The_Function? in a cell to get the description of Name_Of_The_Function and Name_Of_The_Function?? to get the code of the function. Below you find the functions.

LibreNorm(a)
LibreConj(a,b)
revertirF(a)
revertirT(b)

LibreTrenza(libre,trenza)

conjtrenza(lista,n)

relstrenzaconj(lista,n)

relstrenza(lista,n)

invertirlista(lista)

cambio_rel(lista,nuevo,elim)

cambio(grupo,elim_nuevo,elim,lista=[])

cambio_elim(lista,pal,elim)

eliminar(grupo,generador,lista=[])

CyclicComm(lista)

GtoK(G,elto,K)

GtoTietze(G,elto)

abelianizar(tz,m)

MatrizAbel(grupo)

CambioVarSmith(matriz,smith,tt)

caracter(x,R,cambio)

unidades(f,R)

unidadeslista(lista,R)

unidadesmatriz(matriz,R,O='R')

Hay_unidades(A,R,S,dividir=False)

reducir_matriz(A,R,S,dividir=False)

grafoplumbing(grafo,selfint,generos=None,flechas=None)

AnilloMatriz(grupo)

pseudo_coxeter(lista)

gapperm(listatupla)

generadores_trenzas(lista,hilos)

GapConvert(objeto)

HomomorphismGroups(G1,G2,im)

HomImageEl(hom,el)

HomImageSub(hom,lista)

Centralizador(ggap,elgap)

Normalizador(ggap,subgap)

Hurwitz0(B,gensS,br1,br2,GS,lista)

HurwitzTrenza(trenza,lista,grupo)

Hurwitz1(B,gensS,generadores,br1,br2,GS,lista)

HurwitzOrbita(B,gensS,generadores,br1,GS,lista)

burauperm(nm,m,hilos,field=False)

compruebaorbita(orbita,lista)

orbitapermutacion(orbita,grupo,generadores)
