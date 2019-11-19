"""
Created by José Hernán Ortiz Ocampo for the course Optics Instrumentation, 
university EAFIT, 2019-2

Module with the functions corresponding to SistemaEncriptacion_simplus_main

Has the functions:
    
    FT
    iFT
    padding
    encryptation
    desencryptation
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def FT(u):
    """
    Calculates de 2d fft of u, with the corresponding shifts
    """
    U = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u)))
    return U

def iFT(u):
    """
    Calculates de 2d ifft of u, with the corresponding shifts
    """
    U = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(u)))
    return U


def padding(image,n,m,M):
    """
    returns numpy MxM matrix with image with size 2*w (square) in the
    middle of the matrix, and 0s elsewhere
    
    Inputs:
        ui = original matrix, size MxM 
        image = path to image to joint
        L = total size of the image
        w = half width of the image to be assembled
    """
    im = Image.open(image).convert('L')
    u = np.zeros((M,M), dtype=complex)
    im_array = np.array(im)
    im_array = im_array/255
    u[int(M/2-m/2):int(M/2-m/2) + m , int(M/2-n/2):int(M/2-n/2) + n ] = im_array
    return u 
 
    
#%% Encrypting function
def encryptation(image,wi, f1, f2, lambda_ , maskpath, px, doublemask = False) :
    """
    Encryption function, 
    
    Inputs:
        image : path to the image to be encrypted
        wi : size of the image to be encrypted ( if rectangular, use the biggest side)
        f1 : focal 1 of the 4f encrypting system
        f2:  focal 2 of the 4f encrypting system
        lambda_ : wavelenght used in the system
        maskpath : path to the mask image
        px : pixel size of the mask 
        doublemask : activate for double random mask encryption
        
    Outputs:
        
        uo = encrypted wave
        um_phase = random phase introduced
        xo = exit coordinates
    
    """
    mask = Image.open(maskpath).convert('L')
    n,m = mask.size
    mask = np.array(mask)
    M = max(n,m)
    mask_array = np.zeros((M, M))
    mask_array[int(M/2-m/2):int(M/2-m/2) + m , int(M/2-n/2):int(M/2-n/2) + n ] = mask
    mask_size = px*M
    xm= np.array([-mask_size/2, mask_size/2])    
    fx = xm/(lambda_*f1)     
    im = Image.open(image)
    n,m = im.size

    xi = np.array([-wi/2, wi/2])
    xo = f1/f2 * xi #the minus sign goes to the inverted image
    
    if doublemask:
        ui = padding(image,n,m,M) * np.exp(1j*np.random.rand(M,M) * 2*np.pi)
    else:
        ui = padding(image,n,m,M)
        
    u3 = FT(ui) 
    um_phase = mask/255 * 2*np.pi
    um = np.exp(1j*um_phase) # only phase mask
    u4 = u3*um # mask imposition to the wave
    uo = iFT(u4) 
    uo = uo[::-1, ::-1]
    

    # Source wave
    plt.figure()
    uiplot = np.abs(ui[int(M/2-m/2):int(M/2-m/2) + m , int(M/2-n/2):int(M/2-n/2) + n ]**2)
    plt.imshow(uiplot,extent = [xi[0], xi[-1], xi[0], xi[-1]],  cmap = 'gray')
    plt.title('Input image')
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')
    
    #randon phase mask
    plt.figure()
    plt.imshow(np.angle(um),extent = [fx[0], fx[-1], fx[0], fx[-1]], cmap = 'gray')
    plt.title(' Random phase mask' )
    plt.xlabel('fx(cyc/mm)')
    plt.ylabel('fy(cyc/mm)')
    
    plt.figure()
    plt.imshow(np.angle(um),extent = [xm[0], xm[-1], xm[0], xm[-1]], cmap = 'gray')
    plt.title(' Random phase mask' )
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')
    # Output plane
    plt.figure()
    uoplot = np.abs(uo[int(M/2-m/2):int(M/2-m/2) + m , int(M/2-n/2):int(M/2-n/2) + n ]**2)
    plt.imshow(uoplot,extent = [xo[0], xo[-1], xo[0], xo[-1]],  cmap = 'gray')
    plt.title('Encrypted image')
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')

    return uo, um_phase, xo

#%%Decryption  function
def desencryptation(uo,um_phase,n,m, xo, f3, f4):
    """
    Encryption function, 
    
    Inputs:
        uo = encrypted wave
        um_phase = encryption mask used,
        n,m = number of pixels of the original image,
        xo = coordinate position of the encrypted image
        f3, f4 = focal of the decrytion 4f system
        
    Outputs:
        
        uo_deconv = decryted image
    
    """
    
    M,M = uo.shape
    Fuo = iFT(uo)
    Fuo_deconv = Fuo * np.exp(-1j*um_phase)
    uo_deconv = iFT(Fuo_deconv)
    x_deconv = f3/f4 * xo 
    
    plt.figure()
    uo_deconvplot = np.abs(uo_deconv[int(M/2-m/2):int(M/2-m/2) + m , int(M/2-n/2):int(M/2-n/2) + n ]**2)
    plt.imshow(uo_deconvplot,extent = [x_deconv[0], x_deconv[-1], x_deconv[0], x_deconv[-1]],  cmap = 'gray')
    plt.title('Decrypted image')
    plt.xlabel('x(mm)')
    plt.ylabel('y(mm)')
    return uo_deconv
