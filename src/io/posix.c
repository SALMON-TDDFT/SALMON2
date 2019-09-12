/*
 *  Copyright 2019 SALMON developers
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "config.h"


#if defined(SYSTEM_HAS_POSIX) \
    && defined(SYSTEM_HAS_POSIX_STAT) \
    && defined(SYSTEM_HAS_POSIX_ACCESS) \
    && defined(SYSTEM_HAS_POSIX_MKDIR)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h> /* snprintf */

#if defined(SYSTEM_HAS_PATH_MAX_IN_LIMITS_H)
#include <limits.h>
#elif defined(SYSTEM_HAS_PATH_MAX_IN_LINUX_LIMITS_H)
#include <linux/limits.h>
#else
#define PATH_MAX 1024 /* FIXME */
#endif

/*
 * check file exists by POSIX
 *
 * retcode: 0: file not exist
 *          1: file exist
 */
void posix_file_exists(char const* filepath, int * retcode) {
  if (access(filepath, F_OK) == 0) {
    *retcode = 1;
  } else {
    *retcode = 0;
  }
}

/*
 * check directory exists by POSIX
 *
 * retcode: 0: directory not exist
 *          1: directory exist
 */
void posix_directory_exists(char const* dirpath, int * retcode) {
  struct stat buf;
  if (stat(dirpath, &buf) == 0 && S_ISDIR(buf.st_mode)) {
    *retcode = 1;
  } else {
    *retcode = 0;
  }
}

/*
 * create directory by POSIX
 * it create a directory recursively.
 *
 * retcode:         0: success
 *          otherwise: error
 */
void posix_mkdir(char const* dirpath, int * retcode) {
  /* recursive directory create */
  char path_tmp[PATH_MAX];
  size_t len;
  snprintf(path_tmp, sizeof(path_tmp), "%s", dirpath);
  len = strlen(path_tmp);
  if(path_tmp[len - 1] == '/')
    path_tmp[len - 1] = 0;
  for(char *p = path_tmp + 1; *p; p++) {
    if(*p == '/') {
      *p = 0;
      mkdir(path_tmp, 0755);
      *p = '/';
    }
  }
  *retcode = mkdir(path_tmp, 0755);
}

#endif

#if defined(SYSTEM_HAS_STDIO_REMOVE)

#include <stdio.h>

/*
 * remove file or directory
 *
 * retcode:         0: success
 *          otherwise: error
 */
void stdio_remove(char const* dirpath, int * retcode) {
  *retcode = remove(dirpath);
}

#endif
