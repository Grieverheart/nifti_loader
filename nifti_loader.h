#ifndef NT_NIFTI_LOADER_H
#define NT_NIFTI_LOADER_H

#include <stdint.h>
#include <stddef.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#pragma pack(push, 1)
typedef struct
{
    int     sizeof_hdr;      // Size of the header. Must be 348 (bytes).
    char    data_type[10];   // Not used; compatibility with analyze.
    char    db_name[18];     // Not used; compatibility with analyze.
    int     extents;         // Not used; compatibility with analyze.
    short   session_error;   // Not used; compatibility with analyze.
    char    regular;         // Not used; compatibility with analyze.
    char    dim_info;        // Encoding directions (phase, frequency, slice).
    short   dim[8];          // Data array dimensions.
    float   intent_p1;       // 1st intent parameter.
    float   intent_p2;       // 2nd intent parameter.
    float   intent_p3;       // 3rd intent parameter.
    short   intent_code;     // nifti intent.
    short   datatype;        // Data type.
    short   bitpix;          // Number of bits per voxel.
    short   slice_start;     // First slice index.
    float   pixdim[8];       // Grid spacings (unit per dimension).
    float   vox_offset;      // Offset into a .nii file.
    float   scl_slope;       // Data scaling, slope.
    float   scl_inter;       // Data scaling, offset.
    short   slice_end;       // Last slice index.
    char    slice_code;      // Slice timing order.
    char    xyzt_units;      // Units of pixdim[1..4].
    float   cal_max;         // Maximum display intensity.
    float   cal_min;         // Minimum display intensity.
    float   slice_duration;  // Time for one slice.
    float   toffset;         // Time axis shift.
    int     glmax;           // Not used; compatibility with analyze.
    int     glmin;           // Not used; compatibility with analyze.
    char    descrip[80];     // Any text.
    char    aux_file[24];    // Auxiliary filename.
    short   qform_code;      // Use the quaternion fields.
    short   sform_code;      // Use of the affine fields.
    float   quatern_b;       // Quaternion b parameter.
    float   quatern_c;       // Quaternion c parameter.
    float   quatern_d;       // Quaternion d parameter.
    float   qoffset_x;       // Quaternion x shift.
    float   qoffset_y;       // Quaternion y shift.
    float   qoffset_z;       // Quaternion z shift.
    float   srow_x[4];       // 1st row affine transform
    float   srow_y[4];       // 2nd row affine transform.
    float   srow_z[4];       // 3rd row affine transform.
    char    intent_name[16]; // Name or meaning of the data.
    char    magic[4];        // Magic string.
} NtNiftiHeader;
#pragma pack(pop)

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

typedef struct
{
    NtNiftiHeader* header;
    void*          voxel_data;
    uint16_t       voxel_type;
    uint8_t        dim;
    size_t*        shape;
    float          affine[4*4];
} NtNifti;

extern NtNifti nt_nifti_file_read(const char* filepath);
extern void nt_nifti_free(NtNifti* nifti);

#ifdef __cplusplus
}
#endif

#ifdef NT_NIFTI_LOADER_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#pragma pack(push, 1)
typedef struct
{
    uint8_t  magic[2];
    uint8_t  compression_method;
    uint8_t  file_flags;
    uint32_t timestamp;
    uint8_t  compression_flags;
    uint8_t  os_id;
} NtiGzipHeader;
#pragma pack(pop)

enum NtiGzipFileFlagsEnum
{
    GZ_FTEXT    = 0x01, //If set the uncompressed data needs to be treated as text instead of binary data.
                        //This flag hints end-of-line conversion for cross-platform text files but does not enforce it.
    GZ_FHCRC    = 0x02, //The file contains a header checksum (CRC-16)
    GZ_FEXTRA   = 0x04, //The file contains extra fields
    GZ_FNAME    = 0x08, //The file contains an original file name string
    GZ_FCOMMENT = 0x10, //The file contains comment
};


#ifndef NTIW_ASSERT
#include <assert.h>
#define NTIW_ASSERT(x) assert(x)
#endif

#if 0
static const char* NtiGzipFileFlags[] =
{
    "GZ_FTEXT",
    "GZ_FHCRC",
    "GZ_FEXTRA",
    "GZ_FNAME",
    "GZ_FCOMMENT"
};

static const char* NtiGzipOS[] =
{
    "FAT filesystem (MS-DOS, OS/2, NT/Win32)",
    "Amiga",
    "VMS (or OpenVMS)",
    "Unix",
    "VM/CMS",
    "Atari TOS",
    "HPFS filesystem (OS/2, NT)",
    "Macintosh",
    "Z-System",
    "CP/M",
    "TOPS-20",
    "NTFS filesystem (NT)",
    "QDOS",
    "Acorn RISCOS",
    "Unknown",
    "Unknown",
    "BE OS",
    "Unknown",
    "OS/400",
    "Mac OS"
};
static void nti__print_gzip_header(NtiGzipHeader header)
{
    printf("magic: 0x%x 0x%x\n", header.magic[0], header.magic[1]);
    printf("timestamp: %d\n", header.timestamp);
    printf("file flags: 0x%x\n", header.file_flags);
    printf("Compression method: 0x%x\n", header.compression_method);
    printf("Compression flags: 0x%x\n", header.compression_flags);
    for(size_t i = 0; i < 5; ++i)
    {
        if(header.file_flags & (1 << i))
            printf("    File flag %s active\n", NtiGzipFileFlags[i]);
    }
    printf("OS id: 0x%x\n", header.os_id);
    printf("    %s\n", header.os_id < 20? NtiGzipOS[header.os_id]: "Unknown");
}
#endif

typedef struct
{
    size_t byte_id;
    uint8_t bit_id;
    const uint8_t* data;
    size_t size;
} NtiDeflateState;

inline static void nti__deflate_eat_bits(NtiDeflateState* state, uint8_t nbits)
{
    // @todo: Improve!
    state->bit_id += nbits;
    while(state->bit_id >= 8)
    {
        state->bit_id -= 8;
        ++state->byte_id;
    }
}

inline static void nti__deflate_peek_bytes(const NtiDeflateState* state, uint8_t* bytes, uint8_t nbytes)
{
    if(state->bit_id > 0)
    {
        for(uint8_t bi = 0; bi < nbytes; ++bi)
        {
            bytes[bi] = state->data[state->byte_id + bi] >> state->bit_id;
            size_t next_byte_pos = state->byte_id + bi + 1;
            if(next_byte_pos < state->size)
            {
                uint8_t extension = ((uint16_t) state->data[next_byte_pos] << 8) >> state->bit_id;
                bytes[bi] |= extension;
            }
            else break;
        }
    }
    else
    {
        for(uint8_t bi = 0; bi < nbytes; ++bi)
        {
            if(state->byte_id + bi == state->size) break;
            bytes[bi] = state->data[state->byte_id + bi];
        }
    }
}

inline static void nti__deflate_skip_byte(NtiDeflateState* state)
{
    if(state->bit_id > 0)
    {
        state->bit_id = 0;
        ++state->byte_id;
    }
}

inline static uint8_t nti__deflate_read_nbits8(NtiDeflateState* state, uint8_t nbits)
{
    if(!nbits) return 0;

    uint8_t bits = state->data[state->byte_id] >> state->bit_id;
    if(state->bit_id > (8 - nbits))
    {
        ++state->byte_id;
        int num_bits = state->bit_id - (8 - nbits);
        int bit_mask = (1 << num_bits) - 1;
        bits |= (state->data[state->byte_id] & bit_mask) << (nbits - num_bits);
    }
    state->bit_id = (state->bit_id + nbits) % 8;
    if(state->bit_id == 0) ++state->byte_id;

    return bits & ((1 << nbits) - 1);
}

inline static uint16_t nti__deflate_read_nbits16(NtiDeflateState* state, uint8_t nbits)
{
    if(nbits <= 8)
        return (uint16_t) nti__deflate_read_nbits8(state, nbits);

    uint16_t bits_a = (uint16_t) nti__deflate_read_nbits8(state, 8);
    uint16_t bits_b = (uint16_t) nti__deflate_read_nbits8(state, nbits - 8);
    return bits_a | (bits_b << 8);
}

#define __NTI_NB 7
#define BITMASK(x) ((1ULL << (x)) - 1)

struct NtiBinaryTree_t;
typedef struct NtiBinaryTree_t
{
    uint16_t* data;
    uint8_t* lens;
    size_t size_data;
    size_t far_size;
} NtiBinaryTree;

// @todo: Instead of far trees, create subtables. The first __NTI_NB bits of
// the code when indexed in the root table, will point to the subtable.
// Building the subtables so that they can be accessed efficiently is a bit
// tricky; look into libdeflate, specifically in build_decode_table.
inline static void nti__bt_init(NtiBinaryTree* bt)
{
    bt->size_data = 1ULL << __NTI_NB;

    size_t rest_bits = 15 - __NTI_NB;
    bt->far_size = (1ULL << (rest_bits + 1)) - 1;
    size_t num_far_trees = 1ULL << __NTI_NB;
    size_t total_size = (bt->size_data + num_far_trees * bt->far_size);

    bt->data = (uint16_t*) malloc(total_size * sizeof(*bt->data));
    bt->lens = (uint8_t*) malloc(bt->size_data * sizeof(*bt->lens));
    for(size_t i = 0; i < total_size; ++i)
        bt->data[i] = (uint16_t)(-1);
}

inline static void nti__bt_reset(NtiBinaryTree* bt)
{
    size_t num_far_trees = 1ULL << __NTI_NB;
    size_t total_size = (bt->size_data + num_far_trees * bt->far_size);
    for(size_t i = 0; i < total_size; ++i)
        bt->data[i] = (uint16_t)(-1);
}

inline static void nti__bt_free(NtiBinaryTree* bt)
{
    free(bt->data);
    free(bt->lens);
}

inline static void nti__bt_push(NtiBinaryTree* bt, uint16_t data, uint16_t code, uint8_t code_length)
{
    if(code_length > __NTI_NB)
    {
        // Push using the reverse code.
        size_t pos = 0;
        uint16_t mask = 1 << (code_length - 1);
        uint16_t idx = 0;
        for(size_t i = 0; i < __NTI_NB; ++i)
        {
            if(code & mask)
                idx |= (1ULL << i);
            code <<= 1;
        }

        // Push using the reverse code.
        uint16_t* far_data = bt->data + bt->size_data + idx * bt->far_size;
        pos = 0;
        for(size_t i = 0; i < (size_t)(code_length - __NTI_NB); ++i)
        {
            pos = 2 * pos + 1 + ((code & mask) > 0);
            far_data[pos] = (uint16_t)(-2);
            code <<= 1;
        }
        far_data[pos] = data;
    }
    else
    {
        // Push using the reverse code.
        uint16_t rcode = 0;
        uint16_t mask = 1ULL << (code_length - 1);
        for(size_t i = 0; i < (size_t)code_length; ++i)
        {
            if(code & mask)
                rcode |= (1 << i);
            code <<= 1;
        }

        //printf("___ 0x%x, %u\n", rcode, code_length);
        mask = (1ULL << code_length) - 1;
        for(uint16_t i = 0; i < bt->size_data; ++i)
        {
            if((i & mask) == rcode)
            {
                bt->data[i] = data;
                bt->lens[i] = code_length;
            }
        }
    }
}

inline static uint16_t nti__bt_match(NtiBinaryTree* bt, uint16_t code, uint8_t* match_length)
{
    uint16_t data = bt->data[code & BITMASK(__NTI_NB)];
    if(data < (uint16_t)(-1))
    {
        *match_length = bt->lens[code & BITMASK(__NTI_NB)];
        return data;
    }

    size_t idx = code & BITMASK(__NTI_NB);
    code = code >> __NTI_NB;
    uint16_t* far_data = bt->data + bt->size_data + idx * bt->far_size;

    size_t pos = 0;
    for(size_t i = 0; i < 15 - __NTI_NB; ++i)
    {
        pos = (pos << 1) + 1 + (code & 1);

        if(far_data[pos] < (uint16_t)(-2))
        {
            *match_length = i+1+__NTI_NB;
            return far_data[pos];
        }
        if(far_data[pos] == (uint16_t)(-1)) return (uint16_t)(-1);
        code >>= 1;
    }


    return (uint16_t)(-1);
}

inline static void nti__deflate_build_code_tree(NtiBinaryTree* bt, uint8_t* code_lengths, size_t num_codes)
{
    size_t max_code = 0;
    uint16_t* bl_count = (uint16_t*)calloc(16+1, sizeof(uint16_t));
    for(size_t code = 0; code < num_codes; ++code)
    {
        uint8_t code_length = code_lengths[code];
        NTIW_ASSERT(code_length < 17);
        if(code_length > 0)
        {
            max_code = code;
            ++bl_count[code_length];
        }
    }

    uint16_t code = 0;
    uint16_t next_code[16+1];
    for(uint8_t bits = 1; bits < 16; ++bits)
    {
        code = (code + bl_count[bits-1]) << 1;
        next_code[bits] = code;
    }
    free(bl_count);

    for(size_t n = 0; n <= max_code; ++n)
    {
        uint8_t len = code_lengths[n];
        if(len > 0)
        {
            nti__bt_push(bt, n, next_code[len], len);
            ++next_code[len];
        }
    }
}

typedef struct
{
    uint8_t* data;
    size_t capacity;
    size_t count;
} NtiBuffer;

inline static void nti__buffer_init(NtiBuffer* buffer, size_t initial_capacity)
{
    buffer->capacity = (initial_capacity > 0)? initial_capacity: 32;
    buffer->data = (uint8_t*) malloc(buffer->capacity);
    buffer->count = 0;
}

inline static void nti__buffer_append(NtiBuffer* buffer, uint8_t data)
{
    if(buffer->count + 1 >= buffer->capacity)
    {
        buffer->capacity *= 2;
        buffer->data = (uint8_t*) realloc(buffer->data, buffer->capacity);
    }
    buffer->data[buffer->count++] = data;
}

inline static void nti__buffer_append_v(NtiBuffer* buffer, const uint8_t* data, size_t count)
{
    if(buffer->count + count >= buffer->capacity)
    {
        buffer->capacity = (buffer->count + count) * 2;
        buffer->data = (uint8_t*) realloc(buffer->data, buffer->capacity);
    }
    memcpy(buffer->data + buffer->count, data, count);
    buffer->count += count;
}

inline static void nti__buffer_compact(NtiBuffer* buffer)
{
    buffer->capacity = buffer->count;
    buffer->data = (uint8_t*) realloc(buffer->data, buffer->capacity);
}

inline static void nti__buffer_maybe_resize(NtiBuffer* buffer, size_t count)
{
    if(buffer->count + count < buffer->capacity) return;

    buffer->capacity = (buffer->count + count) * 2;
    buffer->data = (uint8_t*) realloc(buffer->data, buffer->capacity);
}

inline static uint8_t* nti__inflate(const uint8_t* data, size_t size, size_t* size_out, size_t size_inflated)
{
    // DEFLATE
    //
    // Each block of compressed data begins with 3 header bits
    // containing the following data:
    //
    //    first bit       BFINAL
    //    next 2 bits     BTYPE
    //
    // Note that the header bits do not necessarily begin on a byte
    // boundary, since a block does not necessarily occupy an integral
    // number of bytes.
    //
    // BFINAL is set if and only if this is the last block of the data
    // set.
    //
    // BTYPE specifies how the data are compressed, as follows:
    //
    //    00 - no compression
    //    01 - compressed with fixed Huffman codes
    //    10 - compressed with dynamic Huffman codes
    //    11 - reserved (error)
    //
    // The only difference between the two compressed cases is how the
    // Huffman codes for the literal/length and distance alphabets are
    // defined.
    //
    // In all cases, the decoding algorithm for the actual data is as
    // follows:
    //
    //    do
    //       read block header from input stream.
    //       if stored with no compression
    //          skip any remaining bits in current partially
    //             processed byte
    //          read LEN and NLEN (see next section)
    //          copy LEN bytes of data to output
    //       otherwise
    //          if compressed with dynamic Huffman codes
    //             read representation of code trees (see
    //                subsection below)
    //          loop (until end of block code recognized)
    //             decode literal/length value from input stream
    //             if value < 256
    //                copy value (literal byte) to output stream
    //             otherwise
    //                if value = end of block (256)
    //                   break from loop
    //                otherwise (value = 257..285)
    //                   decode distance from input stream
    //
    //                   move backwards distance bytes in the output
    //                   stream, and copy length bytes from this
    //                   position to the output stream.
    //          end loop
    //    while not last block
    static uint8_t alphabet[] = {
        16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
    };
    static const uint16_t length_code_length_base[29] = { // size base for length codes 257..285
        3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
        35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
    };
    static const uint16_t length_code_extra_bits[29] = { // extra bits for length codes 257..285
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
    };
    static const uint16_t dist_code_dist_base[30] = { // offset base for distance codes 0..29
        1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
        257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
        8193, 12289, 16385, 24577
    };
    static const uint16_t dist_code_extra_bits[30] = { // extra bits for distance codes 0..29
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        7, 7, 8, 8, 9, 9, 10, 10, 11, 11,
        12, 12, 13, 13
    };

    NtiBinaryTree btlit_fixed, btdist_fixed;
    btlit_fixed.data  = NULL;
    btdist_fixed.data = NULL;
    btlit_fixed.lens  = NULL;
    btdist_fixed.lens = NULL;

    uint8_t* data_buffer = (uint8_t*) malloc(size_inflated);
    size_t buffer_position = 0;

    NtiBinaryTree bt, btlit, btdist;
    nti__bt_init(&bt);
    nti__bt_init(&btlit);
    nti__bt_init(&btdist);

    NtiDeflateState state;
    state.byte_id = 0;
    state.bit_id  = 0;
    state.data    = data;
    state.size    = size;
    while(true)
    {
        uint8_t header = nti__deflate_read_nbits8(&state, 3);
        uint8_t bfinal = header & 1;
        uint8_t btype  = (header >> 1) & 3;
        if(btype == 0)
        {
            nti__deflate_skip_byte(&state);

            const uint8_t* bytes = state.data + state.byte_id;
            uint16_t len  = *(uint16_t*) bytes;
            uint16_t nlen = *(uint16_t*) (bytes + 2);
            if(len)
            {
                memcpy(data_buffer + buffer_position, bytes + 4, len);
                buffer_position += len;
            }
            state.byte_id += 4 + len;
            NTIW_ASSERT(nlen == (len ^ (uint16_t)(-1)));
        }
        else
        {
            NTIW_ASSERT(btype != 3);

            if(btype == 2)
            {
                uint8_t hlit  = nti__deflate_read_nbits8(&state, 5);
                uint8_t hdist = nti__deflate_read_nbits8(&state, 5);
                uint8_t hclen = nti__deflate_read_nbits8(&state, 4);

                // Read code lengths for code length alphabet.
                {
                    uint8_t num_codes = 19;
                    size_t num_code_lengths = hclen + 4;

                    // Read code lengths.
                    uint8_t* code_lengths = (uint8_t*) calloc(num_codes, 1);
                    {
                        for(uint8_t i = 0; i < num_code_lengths; ++i)
                        {
                            uint8_t code = alphabet[i];
                            uint8_t code_length = nti__deflate_read_nbits8(&state, 3);
                            code_lengths[code] = code_length;
                        }
                    }

                    // Build code length Huffman tree.
                    nti__bt_reset(&bt);
                    nti__deflate_build_code_tree(&bt, code_lengths, num_codes);
                    free(code_lengths);
                }

                // Decode code lengths
                uint8_t* alphabet_code_lengths = (uint8_t*) malloc(hlit + hdist + 258u);
                size_t num_code_lengths = 0;
                {
                    uint16_t previous_letter = (uint16_t)(-1);
                    for(size_t i = 0; i < (hlit + hdist + 258u); ++i)
                    {
                        uint8_t bytes[2];
                        nti__deflate_peek_bytes(&state, bytes, 2);
                        uint8_t match_length = 0;
                        uint16_t letter = nti__bt_match(&bt, *(uint16_t*)bytes, &match_length);
                        NTIW_ASSERT(letter < (uint16_t)(-2));
                        nti__deflate_eat_bits(&state, match_length);

                        if(letter >= 16)
                        {
                            uint8_t repeat;
                            if(letter == 16)
                            {
                                repeat = nti__deflate_read_nbits8(&state, 2) + 3;
                                letter = previous_letter;
                            }
                            else if(letter == 17)
                            {
                                repeat = nti__deflate_read_nbits8(&state, 3) + 3;
                                letter = 0;
                            }
                            else
                            {
                                repeat = nti__deflate_read_nbits8(&state, 7) + 11;
                                letter = 0;
                            }

                            for(size_t n = 0; n < repeat; ++n)
                            {
                                alphabet_code_lengths[num_code_lengths] = (uint8_t) letter;
                                ++num_code_lengths;
                            }
                            i += repeat - 1;
                        }
                        else
                        {
                            alphabet_code_lengths[num_code_lengths] = letter;
                            ++num_code_lengths;
                        }

                        previous_letter = letter;
                    }
                }

                nti__bt_reset(&btlit);
                nti__bt_reset(&btdist);
                nti__deflate_build_code_tree(&btlit, alphabet_code_lengths, hlit+257u);
                nti__deflate_build_code_tree(&btdist, alphabet_code_lengths+hlit+257u, hdist+1u);
                free(alphabet_code_lengths);
            }
            else
            {
                if(!btlit_fixed.data)
                {
                    uint8_t alphabet_code_lengths[288+32];
                    memset(alphabet_code_lengths, 8, 144);
                    memset(alphabet_code_lengths+144, 9, 112);
                    memset(alphabet_code_lengths+256, 7, 24);
                    memset(alphabet_code_lengths+280, 8, 8);
                    memset(alphabet_code_lengths+288, 5, 32);
                    nti__bt_init(&btlit_fixed);
                    nti__bt_init(&btdist_fixed);
                    nti__deflate_build_code_tree(&btlit_fixed, alphabet_code_lengths, 288u);
                    nti__deflate_build_code_tree(&btdist_fixed, alphabet_code_lengths+288u, 32u);
                }

                btlit  = btlit_fixed;
                btdist = btdist_fixed;
            }

            while(true)
            {
                uint8_t bytes[2];
                nti__deflate_peek_bytes(&state, bytes, 2);
                uint8_t match_length = 0;
                uint16_t letter = nti__bt_match(&btlit, *(uint16_t*)bytes, &match_length);
                NTIW_ASSERT(letter < (uint16_t)(-2));
                nti__deflate_eat_bits(&state, match_length);

                if(letter == 256u)
                    break;
                else if(letter < 256u)
                    data_buffer[buffer_position++] = letter;
                else
                {
                    letter -= 257u;
                    uint16_t base = length_code_length_base[letter];
                    uint16_t num_extra_bits = length_code_extra_bits[letter];
                    uint16_t length = base + (uint16_t)nti__deflate_read_nbits8(&state, num_extra_bits);

                    nti__deflate_peek_bytes(&state, bytes, 2);
                    uint8_t match_length = 0;
                    uint16_t letter = nti__bt_match(&btdist, *(uint16_t*)bytes, &match_length);
                    NTIW_ASSERT(letter < 30);
                    nti__deflate_eat_bits(&state, match_length);

                    base = dist_code_dist_base[letter];
                    num_extra_bits = dist_code_extra_bits[letter];
                    uint16_t dist = base + nti__deflate_read_nbits16(&state, num_extra_bits);
                    NTIW_ASSERT(buffer_position >= dist);
                    // @note: Slower
                    //size_t pos = data_buffer.count - dist;
                    //for(size_t i = 0; i < length; ++i)
                    //    nti__buffer_append(&data_buffer, data_buffer.data[pos + i % dist]);


                    const uint8_t* data_lookback = data_buffer + buffer_position - dist;
                    if(length < dist)
                    {
                        memcpy(data_buffer + buffer_position, data_lookback, length);
                        buffer_position += length;
                    }
                    else
                    {
                        size_t l = length;
                        while(l >= dist)
                        {
                            memcpy(data_buffer + buffer_position, data_lookback, dist);
                            buffer_position += dist;
                            l -= dist;
                        }
                        if(l > 0)
                        {
                            memcpy(data_buffer + buffer_position, data_lookback, l);
                            buffer_position += l;
                        }
                    }
                }
            }
        }

        if(bfinal) break;
    }

    if(btlit_fixed.data)
    {
        nti__bt_free(&btlit_fixed);
        nti__bt_free(&btdist_fixed);
    }
    nti__bt_free(&btlit);
    nti__bt_free(&btdist);
    nti__bt_free(&bt);

    *size_out = size_inflated;
    return data_buffer;
}

static inline uint32_t nti__crc32(const uint8_t data[], size_t data_length);

static uint8_t* nti__read_gzip(FILE* fp, size_t* data_size)
{
    fseek(fp, 0, SEEK_END);
    size_t file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    NtiGzipHeader header;
    fread(&header, sizeof(NtiGzipHeader), 1, fp);
    if(*(uint16_t*)header.magic != 0x8b1f)
    {
        fseek(fp, 0, SEEK_SET);
        return NULL;
    }

    if(header.file_flags & GZ_FEXTRA)
    {
        NTIW_ASSERT(false);
        // @todo: implement
        //uint16_t xlen;
        //fread(&xlen, 2, 1, fp);
    }


    if(header.file_flags & GZ_FNAME)
    {
        //NtiBuffer buffer;
        //nti__buffer_init(&buffer, 64);
        while(true)
        {
            uint8_t byte = 0;
            fread(&byte, 1, 1, fp);
            //nti__buffer_append(&buffer, byte);
            if(byte == 0) break;
        }
        //nti__buffer_compact(&buffer);
        //free(buffer.data);
    }

    if(header.file_flags & GZ_FCOMMENT)
    {
        //NtiBuffer buffer;
        //nti__buffer_init(&buffer, 64);
        while(true)
        {
            uint8_t byte = 0;
            fread(&byte, 1, 1, fp);
            //nti__buffer_append(&buffer, byte);
            if(byte == 0) break;
        }
        //nti__buffer_compact(&buffer);
        //free(buffer.data);
    }

    if(header.file_flags & GZ_FHCRC)
    {
        // @todo: implement
        NTIW_ASSERT(false);
    }

    //nti__print_gzip_header(header);

    size_t read_bytes = ftell(fp);
    size_t cdata_size = file_size - read_bytes - 8;
    uint8_t* cdata = (uint8_t*) malloc(cdata_size);
    fread(cdata, cdata_size, 1, fp);
    uint32_t crc32, deflated_data_size;
    fread(&crc32, 4, 1, fp);
    fread(&deflated_data_size, 4, 1, fp);

    uint8_t* data = nti__inflate(cdata, cdata_size, data_size, deflated_data_size);
    free(cdata);
    //uint32_t read_crc32 = nti__crc32(data, *data_size);

    //NTIW_ASSERT(read_crc32 == crc32);

    return data;
}

#if 0
static void nti__print_nifti_header(NtNiftiHeader header)
{
    printf("\nNifti Header:\n");
    printf("===========================\n");
    printf("sizeof_hdr: %d\n", header.sizeof_hdr);
    printf("data_type: %s\n", header.data_type);
    printf("db_name: %s\n", header.db_name);
    printf("extents: %d\n", header.extents);
    printf("session_error: %d\n", header.session_error);
    printf("regular: %d\n", header.regular);
    printf("dim_info: %d\n", header.dim_info);
    printf("dim: %d %d %d %d %d %d %d %d\n", header.dim[0], header.dim[1], header.dim[2], header.dim[3], header.dim[4], header.dim[5], header.dim[6], header.dim[7]);
    printf("intent_p1: %f\n", header.intent_p1);
    printf("intent_p2: %f\n", header.intent_p2);
    printf("intent_p3: %f\n", header.intent_p3);
    printf("intent_code: %d\n", header.intent_code);
    printf("datatype: %d\n", header.datatype);
    printf("bitpix: %d\n", header.bitpix);
    printf("slice_start: %d\n", header.slice_start);
    printf("pixdim: %f %f %f %f %f %f %f %f\n", header.pixdim[0], header.pixdim[1], header.pixdim[2], header.pixdim[3], header.pixdim[4], header.pixdim[5], header.pixdim[6], header.pixdim[7]);
    printf("vox_offset: %f\n", header.vox_offset);
    printf("scl_slope: %f\n", header.scl_slope);
    printf("scl_inter: %f\n", header.scl_inter);
    printf("slice_end: %d\n", header.slice_end);
    printf("slice_code: %d\n", header.slice_code);
    printf("xyzt_units: %d\n", header.xyzt_units);
    printf("cal_max: %f\n", header.cal_max);
    printf("cal_min: %f\n", header.cal_min);
    printf("slice_duration: %f\n", header.slice_duration);
    printf("toffset: %f\n", header.toffset);
    printf("glmax: %d\n", header.glmax);
    printf("glmin: %d\n", header.glmin);
    printf("descrip: %s\n", header.descrip);
    printf("aux_file: %s\n", header.aux_file);
    printf("qform_code: %d\n", header.qform_code);
    printf("sform_code: %d\n", header.sform_code);
    printf("quatern_b: %f\n", header.quatern_b);
    printf("quatern_c: %f\n", header.quatern_c);
    printf("quatern_d: %f\n", header.quatern_d);
    printf("qoffset_x: %f\n", header.qoffset_x);
    printf("qoffset_y: %f\n", header.qoffset_y);
    printf("qoffset_z: %f\n", header.qoffset_z);
    printf("srow_x: %f %f %f %f\n", header.srow_x[0], header.srow_x[1], header.srow_x[2], header.srow_x[3]);
    printf("srow_y: %f %f %f %f\n", header.srow_y[0], header.srow_y[1], header.srow_y[2], header.srow_y[3]);
    printf("srow_z: %f %f %f %f\n", header.srow_z[0], header.srow_z[1], header.srow_z[2], header.srow_z[3]);
    printf("intent_name: %s\n", header.intent_name);
    printf("magic: %s\n", header.magic);
}
#endif

void nt_nifti_free(NtNifti* nifti)
{
    free(nifti->header);
    free(nifti->shape);
}

NtNifti nt_nifti_file_read(const char* filepath)
{
    FILE* fp = fopen(filepath, "rb");

    int magic;
    fread(&magic, 1, 4, fp);
    fseek(fp, 0, SEEK_SET);

    size_t data_size;
    uint8_t* data;
    if((magic & 0xFFFF) == 0x8b1f)
    {
        data_size = 0;
        data = nti__read_gzip(fp, &data_size);
    }
    else if(magic == 348)
    {
        fseek(fp, 0, SEEK_END);
        size_t file_size = ftell(fp);
        fseek(fp, 0, SEEK_SET);

        data = (uint8_t*) malloc(file_size);
        fread(data, 1, file_size, fp);
    }
    else
        NTIW_ASSERT(false);

    fclose(fp);

    NtNiftiHeader header = *(NtNiftiHeader*)data;
    //nti__print_nifti_header(header);

    NTIW_ASSERT(header.sizeof_hdr == 348);
    //NTIW_ASSERT(header.qform_code == 1);
    //NTIW_ASSERT(header.sform_code == 1);
    NtNifti nifti;
    nifti.header = (NtNiftiHeader*)data;
    nifti.voxel_type = header.datatype;
    // @todo: Is the offset correct?
    nifti.voxel_data = (void*)(data + (size_t)header.vox_offset);
    nifti.dim = header.dim[0];
    nifti.shape = (size_t*) malloc(nifti.dim * sizeof(*nifti.shape));
    for(size_t i = 0; i < nifti.dim; ++i)
        nifti.shape[i] = header.dim[i+1];

    if(header.sform_code > 0)
    {
        memcpy(nifti.affine, header.srow_x, 4 * sizeof(float));
        memcpy(nifti.affine + 4, header.srow_y, 4 * sizeof(float));
        memcpy(nifti.affine + 8, header.srow_z, 4 * sizeof(float));
    }
    else if(header.qform_code > 0)
    {
        float b = header.quatern_b;
        float c = header.quatern_c;
        float d = header.quatern_d;
        float a = sqrt(1.0 - b*b - c*c - d*d);
        nifti.affine[0]  = a*a+b*b-c*c-d*d;
        nifti.affine[1]  = 2*(b*c-a*d);
        nifti.affine[2]  = 2*(b*d+a*c);
        nifti.affine[3]  = header.qoffset_x;
        nifti.affine[4]  = 2*(b*c+a*d);
        nifti.affine[5]  = a*a+c*c-b*b-d*d;
        nifti.affine[6]  = 2*(c*d-a*b);
        nifti.affine[7]  = header.qoffset_y;
        nifti.affine[8]  = 2*(b*d-a*c);
        nifti.affine[9]  = 2*(c*d+a*b);
        nifti.affine[10] = a*a+d*d-b*b-c*c;
        nifti.affine[11] = header.qoffset_z;

        for(size_t i = 0; i < 3; ++i)
            for(size_t j = 0; j < 3; ++j)
                nifti.affine[4*j+i] *= header.pixdim[1+i];
    }
    else
    {
        // @todo
        NTIW_ASSERT(false);
    }

    nifti.affine[12] = 0.0;
    nifti.affine[13] = 0.0;
    nifti.affine[14] = 0.0;
    nifti.affine[15] = 1.0;

    // @todo
    // * Apply scl_slope and scl_inter.
    // * Use affine fallbacks.
    // * Handle different data types.
    // * Handle endianness by swapping bytes when needed.

    return nifti;
}

static uint32_t NtiCRCTable[] = {
    0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA,
    0x076DC419, 0x706AF48F, 0xE963A535, 0x9E6495A3,
    0x0EDB8832, 0x79DCB8A4, 0xE0D5E91E, 0x97D2D988,
    0x09B64C2B, 0x7EB17CBD, 0xE7B82D07, 0x90BF1D91,
    0x1DB71064, 0x6AB020F2, 0xF3B97148, 0x84BE41DE,
    0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
    0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC,
    0x14015C4F, 0x63066CD9, 0xFA0F3D63, 0x8D080DF5,
    0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172,
    0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B,
    0x35B5A8FA, 0x42B2986C, 0xDBBBC9D6, 0xACBCF940,
    0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59,
    0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116,
    0x21B4F4B5, 0x56B3C423, 0xCFBA9599, 0xB8BDA50F,
    0x2802B89E, 0x5F058808, 0xC60CD9B2, 0xB10BE924,
    0x2F6F7C87, 0x58684C11, 0xC1611DAB, 0xB6662D3D,
    0x76DC4190, 0x01DB7106, 0x98D220BC, 0xEFD5102A,
    0x71B18589, 0x06B6B51F, 0x9FBFE4A5, 0xE8B8D433,
    0x7807C9A2, 0x0F00F934, 0x9609A88E, 0xE10E9818,
    0x7F6A0DBB, 0x086D3D2D, 0x91646C97, 0xE6635C01,
    0x6B6B51F4, 0x1C6C6162, 0x856530D8, 0xF262004E,
    0x6C0695ED, 0x1B01A57B, 0x8208F4C1, 0xF50FC457,
    0x65B0D9C6, 0x12B7E950, 0x8BBEB8EA, 0xFCB9887C,
    0x62DD1DDF, 0x15DA2D49, 0x8CD37CF3, 0xFBD44C65,
    0x4DB26158, 0x3AB551CE, 0xA3BC0074, 0xD4BB30E2,
    0x4ADFA541, 0x3DD895D7, 0xA4D1C46D, 0xD3D6F4FB,
    0x4369E96A, 0x346ED9FC, 0xAD678846, 0xDA60B8D0,
    0x44042D73, 0x33031DE5, 0xAA0A4C5F, 0xDD0D7CC9,
    0x5005713C, 0x270241AA, 0xBE0B1010, 0xC90C2086,
    0x5768B525, 0x206F85B3, 0xB966D409, 0xCE61E49F,
    0x5EDEF90E, 0x29D9C998, 0xB0D09822, 0xC7D7A8B4,
    0x59B33D17, 0x2EB40D81, 0xB7BD5C3B, 0xC0BA6CAD,
    0xEDB88320, 0x9ABFB3B6, 0x03B6E20C, 0x74B1D29A,
    0xEAD54739, 0x9DD277AF, 0x04DB2615, 0x73DC1683,
    0xE3630B12, 0x94643B84, 0x0D6D6A3E, 0x7A6A5AA8,
    0xE40ECF0B, 0x9309FF9D, 0x0A00AE27, 0x7D079EB1,
    0xF00F9344, 0x8708A3D2, 0x1E01F268, 0x6906C2FE,
    0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7,
    0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC,
    0xF9B9DF6F, 0x8EBEEFF9, 0x17B7BE43, 0x60B08ED5,
    0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252,
    0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B,
    0xD80D2BDA, 0xAF0A1B4C, 0x36034AF6, 0x41047A60,
    0xDF60EFC3, 0xA867DF55, 0x316E8EEF, 0x4669BE79,
    0xCB61B38C, 0xBC66831A, 0x256FD2A0, 0x5268E236,
    0xCC0C7795, 0xBB0B4703, 0x220216B9, 0x5505262F,
    0xC5BA3BBE, 0xB2BD0B28, 0x2BB45A92, 0x5CB36A04,
    0xC2D7FFA7, 0xB5D0CF31, 0x2CD99E8B, 0x5BDEAE1D,
    0x9B64C2B0, 0xEC63F226, 0x756AA39C, 0x026D930A,
    0x9C0906A9, 0xEB0E363F, 0x72076785, 0x05005713,
    0x95BF4A82, 0xE2B87A14, 0x7BB12BAE, 0x0CB61B38,
    0x92D28E9B, 0xE5D5BE0D, 0x7CDCEFB7, 0x0BDBDF21,
    0x86D3D2D4, 0xF1D4E242, 0x68DDB3F8, 0x1FDA836E,
    0x81BE16CD, 0xF6B9265B, 0x6FB077E1, 0x18B74777,
    0x88085AE6, 0xFF0F6A70, 0x66063BCA, 0x11010B5C,
    0x8F659EFF, 0xF862AE69, 0x616BFFD3, 0x166CCF45,
    0xA00AE278, 0xD70DD2EE, 0x4E048354, 0x3903B3C2,
    0xA7672661, 0xD06016F7, 0x4969474D, 0x3E6E77DB,
    0xAED16A4A, 0xD9D65ADC, 0x40DF0B66, 0x37D83BF0,
    0xA9BCAE53, 0xDEBB9EC5, 0x47B2CF7F, 0x30B5FFE9,
    0xBDBDF21C, 0xCABAC28A, 0x53B39330, 0x24B4A3A6,
    0xBAD03605, 0xCDD70693, 0x54DE5729, 0x23D967BF,
    0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94,
    0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D
};

uint32_t nti__crc32(const uint8_t data[], size_t data_length) {
    uint32_t crc32 = 0xFFFFFFFFu;

    for (size_t i = 0; i < data_length; i++) {
        const uint32_t lookupIndex = (crc32 ^ data[i]) & 0xff;
        crc32 = (crc32 >> 8) ^ NtiCRCTable[lookupIndex];  // CRCTable is an array of 256 32-bit constants
    }

    // Finalize the CRC-32 value by inverting all the bits
    crc32 ^= 0xFFFFFFFFu;
    return crc32;
}

#endif

#endif
