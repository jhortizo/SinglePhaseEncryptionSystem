"""
Created by José Hernán Ortiz Ocampo for the course Optics Instrumentation, 
university EAFIT, 2019-2

Main program of the SistemaEncriptacion, works only with the functions program 
SistemaEncriptacion_simplus_fcns

Calls the functions for encrytion and decrytion, and calculates de encrypted 
image and checks for the decrypted one.

For necessarry inputs and outputs, go to fcns docstring

"""

from SistemaEncriptacion_simplus_fcns import encryptation, desencryptation
from PIL import Image


# input parameters
image = 'test_images/test1_jh.png'
mask_image = 'test_images/mask_functional.png'
wi = 10
f1 = 100
f2 = 200
f3 = 100
f4 = 200
px = 0.015
lambda_= 633e-6


uo, um_phase , xo = encryptation(image,wi, f1, f2,lambda_,mask_image, px)

im = Image.open(image)
n,m = im.size

uo_deconv = desencryptation(uo,um_phase,n,m,xo,f3,f4)