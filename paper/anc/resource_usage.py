"""Portable peak-memory observation; no mathematical routines."""
import os

def peaks():
    if os.name == 'nt':
        import ctypes
        class P(ctypes.Structure):
            _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(n,ctypes.c_size_t) for n in ('peak','working','qpp','qpc','qnp','qnc','commit','peak_commit')]
        ctypes.windll.kernel32.GetCurrentProcess.restype=ctypes.c_void_p
        ctypes.windll.psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
        p=P(); p.cb=ctypes.sizeof(p)
        if not ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(),ctypes.byref(p),p.cb):
            raise OSError('Cannot measure process memory')
        return p.peak,p.peak_commit
    import resource,sys
    n=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    n=int(n if sys.platform=='darwin' else n*1024)
    return n,n
