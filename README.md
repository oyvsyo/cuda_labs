# CUDA LABS
Particle penetration simulation written on C++ with cuda extension.
>For [google colab](https://colab.research.google.com) usage only

**Thanks to Google Colab platform** GPU cloud computing is highly available now - you can run paralell calculus not only with python, but as you want, because colab notebook session creates virtual machine with unix shell and nvcc is preinstaled there.
Also there is an notebook extension for C++ code in cell (not memmory sharing between cels)
### curand_curand.ipynb
Time comparison between uniform random number generation on CPU and GPU
[Google colab notebook](https://colab.research.google.com/drive/17J1jVTBmZ7ePZfyD4XPGSZ96-fS0V529)
> Random number generation on CUDA took more time cause of thread switch time
 
![alt text](https://imgur.com/xYA0PRQ.png "uniform random generation time")
### GEANT_on_GPU.ipynb
Time comparison between particle penetration simulation on CPU and GPU
[Google colab notebook](https://colab.research.google.com/drive/1tIGtICU3budeBFCy1_s9W2d7ldDqSlP0)
![alt text](https://imgur.com/uvkIopH.png "simulation time")
### LABA2_FInal.ipynb
Volodymyr Klavdiienko's UI and parametrization of particle penetration simulation on GPU
[Google colab notebook](https://colab.research.google.com/drive/15yr_zj2BVLVBq-4YROuBpkZRl3RuSz3s)

## Particle penetration labs
lab4.py and lab9.py - particle penetration simulation with OOP-based approach (user classes such: Particle, Atom, Enviroment).
lab9.py imports lab4.py for routine task.
![alt text](https://imgur.com/5ihIor9.jpg "particles penetration")
