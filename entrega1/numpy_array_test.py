import numpy as np

arr1 = np.array([10,1,2,5,6,2,3,8])
arr2 = np.array([5,5,5,5,5,5,5,5])
ind = np.where(arr1 < arr2)[0]
for i in ind:
    print (i)
# Imprime los indices de cuando un item en el
# arr1 es menor al arr2.
# Se puede hacer lo mismo comparando con un 
# valor fijo.
