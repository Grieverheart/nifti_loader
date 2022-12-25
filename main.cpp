#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "crc32.h"

#pragma pack(1)
struct NiftiHeader
{
    int     sizeof_hdr;      // Size of the header. Must be 348 (bytes).
    char    data_type[10];   // Not used; compatibility with analyze.
    char    db_name[18];     // Not used; compatibility with analyze.
    int     extents;         // Not used; compatibility with analyze.
    short   session_error;   // Not used; compatibility with analyze.
    char    regular;         // Not used; compatibility with analyze.
    char    dim_info;        //  Encoding directions (phase, frequency, slice).
    short   dim[8];          // Data array dimensions.
    float   intent_p1;       // 1st intent parameter.
    float   intent_p2;       // 2nd intent parameter.
    float   intent_p3;       // 3rd intent parameter.
    short   intent_code;     // nifti intent.
    short   datatype;        //  Data type.
    short   bitpix;          //  Number of bits per voxel.
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
};

#pragma pack(1)
struct GzipHeader
{
    uint8_t  magic[2];
    uint8_t  compression_method;
    uint8_t  file_flags;
    uint32_t timestamp;
    uint8_t  compression_flags;
    uint8_t  os_id;
};

enum GzipFileFlags
{
    GZ_FTEXT    = 0x01, //If set the uncompressed data needs to be treated as text instead of binary data.
                        //This flag hints end-of-line conversion for cross-platform text files but does not enforce it.
    GZ_FHCRC    = 0x02, //The file contains a header checksum (CRC-16)
    GZ_FEXTRA   = 0x04, //The file contains extra fields
    GZ_FNAME    = 0x08, //The file contains an original file name string
    GZ_FCOMMENT = 0x10, //The file contains comment
};

const char* GzipFileFlags[]
{
    "GZ_FTEXT",
    "GZ_FHCRC",
    "GZ_FEXTRA",
    "GZ_FNAME",
    "GZ_FCOMMENT"
};

const char* GzipOS[]
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

void print_gzip_header(GzipHeader header)
{
    printf("magic: 0x%x 0x%x\n", header.magic[0], header.magic[1]);
    printf("timestamp: %d\n", header.timestamp);
    printf("file flags: 0x%x\n", header.file_flags);
    printf("Compression method: 0x%x\n", header.compression_method);
    printf("Compression flags: 0x%x\n", header.compression_flags);
    for(size_t i = 0; i < 5; ++i)
    {
        if(header.file_flags & (1 << i))
            printf("    File flag %s active\n", GzipFileFlags[i]);
    }
    printf("OS id: 0x%x\n", header.os_id);
    printf("    %s\n", header.os_id < 20? GzipOS[header.os_id]: "Unknown");
}

struct DeflateState
{
    size_t byte_id;
    uint8_t bit_id;
};

// @todo: max length!
void deflate_eat_bytes(DeflateState* state, uint8_t nbytes)
{
    state->byte_id += nbytes;
}

void deflate_eat_bits(DeflateState* state, uint8_t nbits)
{
    state->bit_id += nbits;
    if(state->bit_id >= 8)
    {
        state->bit_id -= 8;
        ++state->byte_id;
    }
}

// @todo: max length!
void deflate_peak_bytes(const DeflateState* state, const uint8_t* data, uint8_t* bytes, uint8_t nbytes)
{
    bytes[0] = data[state->byte_id];
    if(state->bit_id > 0)
    {
        for(uint8_t i = 0; i < nbytes; ++i)
        {
            bytes[i] = data[state->byte_id + i] >> state->bit_id;
            uint8_t extension = ((uint16_t) data[state->byte_id+i+1] << 8) >> state->bit_id;
            bytes[i] |= extension;
        }
    }
    else memcpy(bytes, data + state->byte_id, nbytes);
}

uint8_t deflate_read_nbits8(DeflateState* state, const uint8_t* data, uint8_t nbits)
{
    uint8_t bits = data[state->byte_id] >> state->bit_id;
    if(state->bit_id > (8 - nbits))
    {
        ++state->byte_id;
        int num_bits = state->bit_id - (8 - nbits);
        int bit_mask = (1 << num_bits) - 1;
        bits |= (data[state->byte_id] & bit_mask) << (nbits - num_bits);
    }
    state->bit_id = (state->bit_id + nbits) % 8;
    if(state->bit_id == 0) ++state->byte_id;

    return bits & ((1 << nbits) - 1);
}

uint16_t deflate_read_nbits16(DeflateState* state, const uint8_t* data, uint8_t nbits)
{
    if(nbits <= 8)
        return (uint16_t) deflate_read_nbits8(state, data, nbits);

    uint16_t bits_a = (uint16_t) deflate_read_nbits8(state, data, 8);
    uint16_t bits_b = (uint16_t) deflate_read_nbits8(state, data, nbits - 8);
    return bits_a | (bits_b << 8);
}

struct BinaryTree
{
    uint8_t* data;
    size_t size_data;
};

void bt_init(BinaryTree* bt)
{
    bt->size_data = 31;
    bt->data = (uint8_t*) malloc(bt->size_data);
    for(size_t i = 0; i < bt->size_data; ++i)
        bt->data[i] = 255;
}

void bt_free(BinaryTree* bt)
{
    free(bt->data);
}

void bt_push(BinaryTree* bt, uint8_t data, uint8_t code, uint8_t code_length)
{
    size_t new_size = (1 << (size_t)(code_length + 1)) - 1;
    if(new_size > bt->size_data)
    {
        bt->data = (uint8_t*) realloc((uint8_t*) bt->data, new_size);
        for(size_t i = bt->size_data; i < new_size; ++i)
            bt->data[i] = 255;
        bt->size_data = new_size;
    }

    //for(size_t i = 0; i < (size_t)(8 - code_length); ++i)
    //    printf(" ");
    //for(size_t i = 0; i < (size_t)code_length; ++i)
    //{
    //    size_t temp = code_length - i - 1;
    //    printf("%d", (code & (1 << temp)) > 0);
    //}
    //printf("\n");

    size_t pos = 0;
    uint8_t mask = 1 << (code_length - 1);
    for(size_t i = 0; i < (size_t)code_length; ++i)
    {
        pos = 2 * pos + 1 + ((code & mask) > 0);
        code <<= 1;
        bt->data[pos] = 254;
    }
    bt->data[pos] = data;
}

uint8_t bt_match(BinaryTree* bt, uint16_t code, uint8_t* match_length)
{
    size_t pos = 0;
    for(size_t i = 0; i < 16; ++i)
    {
        pos = 2 * pos + 1 + (code & 1);
        if(pos >= bt->size_data || bt->data[pos] == 255)
            return 255;
        //printf("__ %lu, %d, %u\n", pos, code & 1, bt->data[pos]);
        code >>= 1;
        if(bt->data[pos] < 254)
        {
            *match_length = i+1;
            return bt->data[pos];
        }
    }

    return 255;
}

void deflate(const uint8_t* data, size_t size)
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
      // 6,  0,  0, 7, 5, 5, 3, 5,  2, 5,  2, 5,  3, 5,  5, 0,  7,
    };

    printf("\n");
    printf("Starting deflate...\n");

    DeflateState state;
    state.byte_id = 0;
    state.bit_id  = 0;
    while(true)
    {
        uint8_t header = deflate_read_nbits8(&state, data, 3);
        uint8_t bfinal = header & 1;
        uint8_t btype  = (header >> 1) & 3;
        if(btype == 0)
            printf("No compression not supported yet.\n");
        else if(btype == 3)
            assert(false);
        else
        {
            if(btype == 2)
            {
                printf("Reading code trees...\n");
                uint8_t hlit  = deflate_read_nbits8(&state, data, 5);
                uint8_t hdist = deflate_read_nbits8(&state, data, 5);
                uint8_t hclen = deflate_read_nbits8(&state, data, 4);
                printf("%u, %u, %u\n", 257 + hlit, 1 + hdist, 4 + hclen);

                printf("Reading code lengths...\n");
                uint8_t* code_lengths = (uint8_t*) calloc(19, 1);
                uint8_t* codes = (uint8_t*) malloc(19);
                uint8_t max_code = 0;
                {
                    uint8_t bl_count[16] = {};
                    for(uint8_t i = 0; i < hclen + 4; ++i)
                    {
                        uint8_t code = alphabet[i];
                        uint8_t code_length = deflate_read_nbits8(&state, data, 3);
                        code_lengths[code] = code_length;
                        if(code_length > 0 && code > max_code)
                            max_code = code;
                        ++bl_count[code_length];
                        printf("%u, ", code_length);
                    }
                    printf("\n");

                    bl_count[0] = 0;
                    uint8_t code = 0;
                    uint8_t next_code[16];
                    for(uint8_t bits = 1; bits < 16; ++bits)
                    {
                        code = (code + bl_count[bits-1]) << 1;
                        next_code[bits] = code;
                    }

                    for(uint8_t n = 0;  n <= max_code; ++n)
                    {
                        uint8_t len = code_lengths[n];
                        if(len != 0)
                        {
                            codes[n] = next_code[len];
                            ++next_code[len];
                        }
                    }
                }

                BinaryTree bt;
                bt_init(&bt);

                for(uint8_t n = 0;  n <= max_code; ++n)
                {
                    uint8_t len = code_lengths[n];
                    if(len != 0)
                    {
                        //printf("%u\n", n);
                        bt_push(&bt, n, codes[n], code_lengths[n]);
                    }
                }

                free(code_lengths);
                free(codes);

                uint8_t previous_letter = 255;
                size_t num_code_lengths = 0;
                uint8_t* alphabet_code_lengths = (uint8_t*) malloc(hlit + hdist + 258u);
                for(size_t i = 0; i < (hlit + hdist + 258u); ++i)
                {
                    uint8_t bytes[2];
                    deflate_peak_bytes(&state, data, bytes, 2);
                    uint8_t match_length = 0;
                    uint8_t letter = bt_match(&bt, *(uint16_t*)bytes, &match_length);
                    assert(letter < 254);
                    deflate_eat_bits(&state, match_length);

                    if(letter == 16)
                    {
                        uint8_t repeat = deflate_read_nbits8(&state, data, 2) + 3;
                        letter = previous_letter;
                        for(size_t n = 0; n < repeat; ++n)
                        {
                            alphabet_code_lengths[num_code_lengths] = letter;
                            //printf("litlen %lu %d\n", i+n, previous_letter);
                            ++num_code_lengths;
                        }
                        i += repeat - 1;
                    }
                    else
                    {
                        //if(letter > 0) printf("litlen %lu %d\n", i, letter);
                        alphabet_code_lengths[num_code_lengths] = letter;
                        ++num_code_lengths;
                    }

                    previous_letter = letter;
                }
                //printf("--- %lu, %u\n", num_code_lengths, hlit + hdist + 258u);
                bt_free(&bt);

                free(alphabet_code_lengths);
                break;
            }
        }

        if(bfinal) break;
    }
}

void read_gzip(FILE* fp)
{
    fseek(fp, 0, SEEK_END);
    size_t file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    printf("File size: %lu bytes\n", file_size);
    GzipHeader header;
    printf("Reading %lu bytes of header.\n", sizeof(GzipHeader));
    fread(&header, sizeof(GzipHeader), 1, fp);
    print_gzip_header(header);

    // @todo: Optional headers

    size_t data_size = file_size - 10 - 8;
    uint8_t* cdata = (uint8_t*) malloc(data_size);
    fread(cdata, data_size, 1, fp);
    uint32_t crc32, deflated_data_size;
    fread(&crc32, 4, 1, fp);
    fread(&deflated_data_size, 4, 1, fp);

    printf("\n");
    printf("Data size: %f Mb\n", (data_size / 1024.0) / 1024.0);
    printf("Uncompressed data size: %f Mb\n", (deflated_data_size / 1024.0) / 1024.0);
    printf("Compression ratio: %fx\n", (float)deflated_data_size / data_size);
    printf("CRC32: 0x%x\n", crc32);

    deflate(cdata, data_size);

    free(cdata);
}

int main(int argc, char* argv[])
{
    FILE* fp = fopen(argv[1], "rb");
    read_gzip(fp);
    fclose(fp);
    return 0;
}
