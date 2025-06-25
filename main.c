#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>   // For round()
#include <limits.h> // For INT16_MAX, INT16_MIN, INT32_MAX, INT32_MIN
#include <stdbool.h>
#include "hack_effects.h"

// --- Fixed-Point Configuration ---
// We'll now use a Q23 fixed-point format for 32-bit processing.
// This means: 1 sign bit, 8 integer bits, 23 fractional bits.
// The value range for a Q23 fixed_point_t (int32_t) is approximately -256.0 to +255.999...
// This provides more headroom for intermediate calculations if values temporarily exceed -1.0 to 1.0.
#define Q_FORMAT 23
#define FIXED_POINT_ONE (1 << Q_FORMAT) // Represents 1.0 in Q23 format

// Define the fixed_point_t type as a signed 32-bit integer for processing
typedef int32_t fixed_point_t;

// --- Fixed-Point Utility Functions ---
// Convert a float to fixed-point (Q23)
fixed_point_t float_to_q23(float f) {
    // Clamp the float value to prevent overflow if it's too large for the Q23 range.
    // The theoretical max for Q23 (int32_t) with 8 integer bits is (2^8 - 1) / (2^0) = 255.0
    // and min is -256.0. We'll clamp to a slightly safer range to avoid edge cases.
    if (f >= 256.0f) f = 255.999999f;
    if (f < -256.0f) f = -256.0f;
    return (fixed_point_t)round(f * FIXED_POINT_ONE);
}

// Convert fixed-point (Q23) to float
float q23_to_float(fixed_point_t q) {
    return (float)q / FIXED_POINT_ONE;
}

// Fixed-point multiplication (Q23 * Q23 -> Q23)
fixed_point_t q23_mul(fixed_point_t a, fixed_point_t b) {
    // Multiply two 32-bit numbers, result can be up to 64-bit.
    int64_t result_64 = (int64_t)a * b;
    // Right-shift by Q_FORMAT to bring it back to Q23.
    // This implicitly truncates fractional bits beyond Q_FORMAT.
    fixed_point_t final_result = (fixed_point_t)(result_64 >> Q_FORMAT);

    // Optional: Add saturation/clipping if the result is expected to overflow int32_t
    // For normal Q23*Q23 where operands are within Q23, the result should also fit if
    // the integer part of Q23 has enough bits. Here, we have 8 integer bits, so if
    // two numbers like 100.0 * 100.0 (10000.0) are multiplied, it will overflow
    // the 8 integer bits. Consider adding a saturation here if inputs can be large.
    // For typical audio (normalized -1.0 to 1.0), this is usually fine.
    return final_result;
}

// Fixed-point addition (Q23 + Q23 -> Q23) with saturation
fixed_point_t q23_add(fixed_point_t a, fixed_point_t b) {
    // Use a 64-bit temporary variable to detect overflow before casting back to 32-bit.
    int64_t sum = (int64_t)a + b;

    // Saturate the sum to the valid range of fixed_point_t (int32_t)
    if (sum > INT32_MAX) {
        return INT32_MAX;
    } else if (sum < INT32_MIN) {
        return INT32_MIN;
    }
    return (fixed_point_t)sum;
}

// --- WAV File Header Structure ---
// Standard WAV file header for PCM data
typedef struct {
    // RIFF Chunk
    char     chunkID[4];       // "RIFF"
    uint32_t chunkSize;        // Size of the rest of the file in bytes (fileSize - 8)
    char     format[4];        // "WAVE"

    // FMT Subchunk
    char     subchunk1ID[4];   // "fmt "
    uint32_t subchunk1Size;    // 16 for PCM
    uint16_t audioFormat;      // 1 for PCM
    uint16_t numChannels;      // 1 for Mono, 2 for Stereo
    uint32_t sampleRate;       // 8000, 44100, etc.
    uint32_t byteRate;         // sampleRate * numChannels * bitsPerSample / 8
    uint16_t blockAlign;       // numChannels * bitsPerSample / 8
    uint16_t bitsPerSample;    // 8, 16, 24, 32

    // DATA Subchunk
    char     subchunk2ID[4];   // "data"
    uint32_t subchunk2Size;    // numSamples * numChannels * bitsPerSample / 8
} WavHeader;

// --- Utility function to print header info ---
void print_wav_header_info(const WavHeader *header) {
    printf("WAV Header Information:\n");
    printf("  Chunk ID: %.4s\n", header->chunkID);
    printf("  Chunk Size: %u bytes\n", header->chunkSize);
    printf("  Format: %.4s\n", header->format);
    printf("  Subchunk1 ID: %.4s\n", header->subchunk1ID);
    printf("  Subchunk1 Size: %u\n", header->subchunk1Size);
    printf("  Audio Format: %u (1=PCM)\n", header->audioFormat);
    printf("  Channels: %u\n", header->numChannels);
    printf("  Sample Rate: %u Hz\n", header->sampleRate);
    printf("  Byte Rate: %u bytes/sec\n", header->byteRate);
    printf("  Block Align: %u bytes\n", header->blockAlign);
    printf("  Bits Per Sample: %u\n", header->bitsPerSample);
    printf("  Subchunk2 ID: %.4s\n", header->subchunk2ID);
    printf("  Subchunk2 Size: %u bytes (Audio Data Size)\n", header->subchunk2Size);
    printf("  Total Samples (per channel): %u\n", header->subchunk2Size / (header->numChannels * (header->bitsPerSample / 8)));
}

// --- Function to read WAV file ---
// Reads a 16-bit PCM WAV file, converts to fixed_point_t (Q23) array for processing.
// Returns the number of samples read (total samples, not per channel), or 0 on error.
// The `data` pointer will be allocated inside this function, caller must free it.
uint32_t read_wav_file(const char *filename, fixed_point_t **data, WavHeader *header_out) {
    FILE *file = NULL;
    uint32_t num_samples_total = 0; // Total samples (channel 0 + channel 1 + ...)
    int16_t *raw_pcm_data = NULL;   // Temporary buffer for raw 16-bit PCM data
    fixed_point_t *audio_data_q23 = NULL; // Buffer for Q23 fixed-point data

    printf("Attempting to open WAV file: %s\n", filename);
    file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening input WAV file");
        return 0;
    }

    // Read WAV header
    if (fread(header_out, 1, sizeof(WavHeader), file) != sizeof(WavHeader)) {
        fprintf(stderr, "Error reading WAV header from %s\n", filename);
        fclose(file);
        return 0;
    }

    // Validate header
    if (strncmp(header_out->chunkID, "RIFF", 4) != 0 ||
        strncmp(header_out->format, "WAVE", 4) != 0 ||
        strncmp(header_out->subchunk1ID, "fmt ", 4) != 0 ||
        strncmp(header_out->subchunk2ID, "data", 4) != 0 ||
        header_out->audioFormat != 1 || // Not PCM
        header_out->bitsPerSample != 16) { // Not 16-bit PCM
        fprintf(stderr, "Unsupported WAV format: Must be 16-bit PCM WAV.\n");
        fclose(file);
        return 0;
    }

    print_wav_header_info(header_out);

    // Calculate total number of samples (all channels combined) for 16-bit data
    num_samples_total = header_out->subchunk2Size / (header_out->bitsPerSample / 8);
    printf("Allocating memory for %u raw 16-bit samples...\n", num_samples_total);

    raw_pcm_data = (int16_t *)malloc(num_samples_total * sizeof(int16_t));
    if (!raw_pcm_data) {
        perror("Error allocating memory for raw PCM data");
        fclose(file);
        return 0;
    }

    // Read raw 16-bit PCM audio data
    if (fread(raw_pcm_data, sizeof(int16_t), num_samples_total, file) != num_samples_total) {
        fprintf(stderr, "Error reading raw audio data from %s\n", filename);
        free(raw_pcm_data);
        fclose(file);
        return 0;
    }
    fclose(file); // Close file after reading raw data

    printf("Successfully read %u raw 16-bit samples from %s\n", num_samples_total, filename);

    // Now, convert 16-bit PCM data to Q23 fixed-point format
    printf("Converting raw 16-bit samples to Q23 fixed-point...\n");
    audio_data_q23 = (fixed_point_t *)malloc(num_samples_total * sizeof(fixed_point_t));
    if (!audio_data_q23) {
        perror("Error allocating memory for Q23 audio data");
        free(raw_pcm_data);
        return 0;
    }

    // The conversion from 16-bit signed integer to Q23 fixed-point (32-bit)
    // involves left-shifting the 16-bit value by (Q_FORMAT - 15) bits.
    // This effectively scales the 16-bit fractional representation to the 23-bit fractional.
    for (uint32_t i = 0; i < num_samples_total; ++i) {
        audio_data_q23[i] = ((fixed_point_t)raw_pcm_data[i]) << (Q_FORMAT - 15);
    }

    free(raw_pcm_data); // Free the temporary raw PCM buffer

    *data = audio_data_q23;
    printf("Conversion to Q23 complete.\n");
    return num_samples_total;
}

// --- Function to write WAV file ---
// Writes a fixed_point_t (Q23) array to a 16-bit PCM WAV file.
// Returns 1 on success, 0 on failure.
int write_wav_file(const char *filename, const fixed_point_t *data, uint32_t num_samples_total, const WavHeader *header_in) {
    FILE *file = NULL;
    WavHeader output_header = *header_in; // Copy the input header for modification
    int16_t *output_pcm_data = NULL; // Buffer for 16-bit PCM data to write

    printf("Attempting to open output WAV file: %s\n", filename);
    file = fopen(filename, "wb");
    if (!file) {
        perror("Error opening output WAV file");
        return 0;
    }

    // Update chunk sizes in the header for the output file
    output_header.subchunk2Size = num_samples_total * (output_header.bitsPerSample / 8);
    output_header.chunkSize = 36 + output_header.subchunk2Size; // 36 bytes for header before data

    printf("Writing WAV header to %s...\n", filename);
    print_wav_header_info(&output_header);

    if (fwrite(&output_header, 1, sizeof(WavHeader), file) != sizeof(WavHeader)) {
        fprintf(stderr, "Error writing WAV header to %s\n", filename);
        fclose(file);
        return 0;
    }

    // Convert Q23 fixed-point data back to 16-bit PCM
    printf("Converting Q23 fixed-point samples back to 16-bit PCM...\n");
    output_pcm_data = (int16_t *)malloc(num_samples_total * sizeof(int16_t));
    if (!output_pcm_data) {
        perror("Error allocating memory for output PCM data");
        fclose(file);
        return 0;
    }

    // The conversion from Q23 fixed-point (32-bit) to 16-bit signed integer PCM
    // involves right-shifting by (Q_FORMAT - 15) bits and then saturating.
    for (uint32_t i = 0; i < num_samples_total; ++i) {
        // Shift back to the 16-bit range
        int32_t temp_val = data[i] >> (Q_FORMAT - 15);

        // Saturate the value to the int16_t range [-32768, 32767]
        if (temp_val > INT16_MAX) {
            output_pcm_data[i] = INT16_MAX;
        } else if (temp_val < INT16_MIN) {
            output_pcm_data[i] = INT16_MIN;
        } else {
            output_pcm_data[i] = (int16_t)temp_val;
        }
    }

    // Write audio data
    printf("Writing %u samples to %s...\n", num_samples_total, filename);
    if (fwrite(output_pcm_data, sizeof(int16_t), num_samples_total, file) != num_samples_total) {
        fprintf(stderr, "Error writing audio data to %s\n", filename);
        free(output_pcm_data);
        fclose(file);
        return 0;
    }

    free(output_pcm_data); // Free the temporary output PCM buffer
    fclose(file);
    printf("Successfully wrote processed audio to %s\n", filename);
    return 1;
}

// --- Custom Audio Processing Block (Placeholder) ---
// This function takes an array of fixed_point_t audio samples and applies a simple gain.
// Users can modify this function or add new ones and select them via command-line options.
// `num_samples_total` is the total number of samples across all channels.

void custom_processing_block(fixed_point_t *audio_data, uint32_t num_samples_total, const WavHeader *header) {
    printf("Applying custom processing (simple gain) to %u samples...\n", num_samples_total);

    /* Enable the hack that needs to be tested */
    hack_vinyl_crackle_effect(audio_data, audio_data, num_samples_total);
    hack_ambient_noise_effect(audio_data, audio_data, num_samples_total);
    printf("Custom processing complete.\n");
}

// --- Main function ---
int main(int argc, char *argv[]) {
    fixed_point_t *audio_data = NULL;
    WavHeader input_header;
    uint32_t num_samples_read = 0;

    printf("--- Fixed-Point Audio Processing Utility (int32_t processing) ---\n");

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_wav_file> <output_wav_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_filename = argv[1];
    const char *output_filename = argv[2];

    // 1. Read input WAV file and convert to Q23 fixed-point
    num_samples_read = read_wav_file(input_filename, &audio_data, &input_header);
    if (num_samples_read == 0) {
        fprintf(stderr, "Failed to read input WAV file. Exiting.\n");
        return EXIT_FAILURE;
    }

    // 2. Apply custom processing using Q23 fixed-point arithmetic
    custom_processing_block(audio_data, num_samples_read, &input_header);

    // 3. Write processed audio to output WAV file, converting back to 16-bit PCM
    if (!write_wav_file(output_filename, audio_data, num_samples_read, &input_header)) {
        fprintf(stderr, "Failed to write output WAV file. Exiting.\n");
        free(audio_data); // Free memory even on write failure
        return EXIT_FAILURE;
    }

    // Clean up
    free(audio_data);
    printf("Processing complete. Output saved to %s\n", output_filename);

    return EXIT_SUCCESS;
}
