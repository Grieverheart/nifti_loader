# WIP: NiftiLoader

NiftiLoader is a work in progress (WIP) single header library for loading nifti files. Currently, only gzipped Nifti1 files are supported. Furthermore, the gzip implementation is also not complete.

This library will be expanded as new types of nifti files are encountered while working on other projects that make use of this library.

## Usage

Include the header file `nifti_loader.h` and `#define NT_NIFTI_LOADER_IMPLEMENTATION` in only one translation unit of your choosing.

There are just two functions defined:
``` C
NtNifti nt_nifti_file_read(char* filepath);
void nt_nifti_free(NtNifti* nifti);
```
where
```C
typedef struct
{
    void*       voxel_data;
    uint16_t    voxel_type;
    uint8_t     dim;
    size_t*     shape;
    float       affine[4*4];
    NtNiftiHeader header;
} NtNifti;
```
The `voxel_type` is taken directly from the Nifti1 header and is equal to:
```C
enum NtNiftiTypeEnum
{
    NT_NIFTI_TYPE_UNKNOWN    = 0,
    NT_NIFTI_TYPE_BINARY     = 1,
    NT_NIFTI_TYPE_UINT8      = 2,
    NT_NIFTI_TYPE_INT16      = 4,
    NT_NIFTI_TYPE_INT32      = 8,
    NT_NIFTI_TYPE_FLOAT32    = 16,
    NT_NIFTI_TYPE_COMPLEX64  = 32,
    NT_NIFTI_TYPE_FLOAT64    = 64,
    NT_NIFTI_TYPE_RGB24      = 128,
    NT_NIFTI_TYPE_ALL        = 255,
    NT_NIFTI_TYPE_INT8       = 256,
    NT_NIFTI_TYPE_UINT16     = 512,
    NT_NIFTI_TYPE_UINT32     = 768,
    NT_NIFTI_TYPE_INT64      = 1024,
    NT_NIFTI_TYPE_UINT64     = 1280,
    NT_NIFTI_TYPE_FLOAT128   = 1512,
    NT_NIFTI_TYPE_COMPLEX128 = 1792,
    NT_NIFTI_TYPE_COMPLEX256 = 2048,
    NT_NIFTI_TYPE_RGBA       = 2304
};
```
The goal of the library is to put as much information in the top level of the `NtNifti` struct, which means that the affine transform is calculated based on header information, and the voxel data are also scaled based on the header `scl_slope` and `scl_inter` fields.
