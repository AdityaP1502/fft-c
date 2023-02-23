#ifndef DLL_EXPORT_API_H
  #define DLL_EXPORT_API_H

  #define FFT_EXPORTS

  #ifdef _WIN32

    /* You should define ADD_EXPORTS *only* when building the DLL. */
    #ifdef FFT_EXPORTS
      #define FFTLIBRARY_API __declspec(dllexport)
    #else
      #define FFTLIBRARY_API __declspec(dllimport)
    #endif

    /* Define calling convention in one place, for convenience. */
    #define FFTLIBRARY_CALL __cdecl

  #else /* _WIN32 not defined. */

    /* Define with no value on non-Windows OSes. */
    #define FFTLIBRARY_API
    #define FFTLIBRARY_CALL
    
  #endif
#endif