# Miscellanous-SageFunctions
Some auxiliar functions in **Sagemath**
It contains right now a file named **MonodromyGroupFunctions.sage** with functions mostly related with braid monodromy. Functions are documented. Last update: 11-06-2016.

I use the following code to use these functions in a worksheet or a sage session:

```python
import os
DIR='path_where_the_file_is_universally_accessible/'  
SAGE_LOCATION='path_where_sagemath_is_installed/'  
if 'MonodromyGroupFunctions.sage.py' in os.listdir(DIR):  
    NO=os.system('rm '+DIR+'MonodromyGroupFunctions.sage.py')  
NO=os.system(SAGE_LOCATION+'local/bin/sage-preparse '+DIR+'MonodromyGroupFunctions.sage')  
os.chmod(DIR+'MonodromyGroupFunctions.sage.py',0o666)  
load(DIR+'MonodromyGroupFunctions.sage.py')
```  
 

In a sage worksheet type**Name_Of_The_Function?** in a cell to get the description of 
**Name_Of_The_Function** and **Name_Of_The_Function??** to get the code of the function. Below you find the functions.

1. LibreNorm(a)
2. LibreConj(a,b)
3. revertirF(a)
4. revertirT(b)
5. LibreTrenza(libre,trenza)
6. conjtrenza(lista,n)
7. relstrenzaconj(lista,n)
8. relstrenza(lista,n)
9. invertirlista(lista)
1. cambio_rel(lista,nuevo,elim)
1. cambio(grupo,elim_nuevo,elim,lista=[])
1. cambio_elim(lista,pal,elim)
1. eliminar(grupo,generador,lista=[])
1. CyclicComm(lista)
1. GtoK(G,elto,K)
1. GtoTietze(G,elto)
1. abelianizar(tz,m)
1. MatrizAbel(grupo)
1. CambioVarSmith(matriz,smith,tt)
1. caracter(x,R,cambio)
1. unidades(f,R)
1. unidadeslista(lista,R)
1. unidadesmatriz(matriz,R,O='R')
1. Hay_unidades(A,R,S,dividir=False)
1. reducir_matriz(A,R,S,dividir=False)
1. grafoplumbing(grafo,selfint,generos=None,flechas=None)
1. AnilloMatriz(grupo)
1. pseudo_coxeter(lista)
1. gapperm(listatupla)
1. generadores_trenzas(lista,hilos)
1. GapConvert(objeto)
1. HomomorphismGroups(G1,G2,im)
1. HomImageEl(hom,el)
1. HomImageSub(hom,lista)
1. Centralizador(ggap,elgap)
1. Normalizador(ggap,subgap)
1. Hurwitz0(B,gensS,br1,br2,GS,lista)
1. HurwitzTrenza(trenza,lista,grupo)
1. Hurwitz1(B,gensS,generadores,br1,br2,GS,lista)
1. HurwitzOrbita(B,gensS,generadores,br1,GS,lista)
1. burauperm(nm,m,hilos,field=False)
1. compruebaorbita(orbita,lista)
1. orbitapermutacion(orbita,grupo,generadores)
