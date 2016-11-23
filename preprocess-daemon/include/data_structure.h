#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

#include <csv_file.h>

/*
 * Structure describing the feature vector
 * frame information such as the presence
 * of an eye-blink and such... 
 */
typedef struct frame_info_s{
	char eye_blink_detected;
	char padding[7];
}frame_info_t;

/*
 * Structure describing the feature vector 
 * it is terminated by a frame_status and
 * the actual feature vector
 */
typedef struct feature_buf_s {
	int nb_features;
	unsigned char *featvect_ptr;
	frame_info_t frame_status; // note this field is not supported by the CSV, yet
} feature_buf_t;




#endif
