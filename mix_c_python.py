import numpy as np
import ctypes
from ctypes import *  # 传参用 ctypes.c_char_p


so = ctypes.cdll.LoadLibrary
# huffman.so: 编译得到的动态库
lib = so("./huffman.so")   
# name_p 传入变量的指针, 方便c++代码中获取"xxx.npy"信号
# .encode("utf-8") 应该是py3需要的编码方式 
name_p = ctypes.c_char_p("xxx.npy".encode("utf-8"))
lib.fun(name_p)  # 调用so动态库中的fun()函数, 且传入参数为"xxx.npy"