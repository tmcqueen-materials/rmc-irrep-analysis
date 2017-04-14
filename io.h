// stdout/stderr IO

void printErrorAndExit(int err, char *strErr, ...) {
	va_list args;
	va_start (args, strErr);
	if (strErr != NULL) vfprintf(stderr, strErr, args);
	va_end (args);
	exit(err);
}

void printWarning(char *strWarn, ...) {
        va_list args;
        va_start (args, strWarn);
        if (strWarn != NULL) vfprintf(stderr, strWarn, args);
        va_end (args);
}

void printStatus(char *strStat, ...) {
        va_list args;
        va_start (args, strStat);
        if (strStat != NULL) vfprintf(stdout, strStat, args);
        va_end (args);
}

void printUsageAndExit(int err) {
	printErrorAndExit(err, USAGE_MESSAGE);
}

// BEGIN MINI FILE PARSING MODULE

#define FTA_FILE 0
#define TMP_FILE 1

int sscanf_str(char *str, char *strBuf, int strBufLen) {
	char tmp[MAX_LINE_LEN] = { 0 };
	char fmt[16] = { 0 };
	if (str == NULL || strBuf == NULL || strBufLen < 2) return 0;
	if (sscanf(str, " %s", tmp) < 1) return 0;
	if (snprintf(fmt, 14, "%%%is", strBufLen-1) < 1) return 0;
	return sscanf(tmp, fmt, strBuf);
}

void closeFile(int f) {
	if (f < 0 || f > MAX_OPEN_FILES-1) { return; }
        if (curFile[f] != NULL) { fclose(curFile[f]); curFile[f] = NULL; }
}

bool isOpen(int f) {
        if (f < 0 || f > MAX_OPEN_FILES-1) { return false; }
	if (curFile[f] != NULL) return true;
	return false;
}

void openFile(int f, char *fileName) {
	closeFile(f);
	if (f < 0 || f > MAX_OPEN_FILES-1) printErrorAndExit(2, "Invalid File Handle!\n");
	if (fileName == NULL || strlen(fileName) < 1) printErrorAndExit(2, "Invalid/Null Filename!\n");
	curFile[f] = fopen(fileName, "rb");
	if (!curFile[f]) printErrorAndExit(3, "Could not open file %s!\n", fileName);
	fprintf(stdout, "Reading %s...\n", fileName);
}

// obeys perl's convention that lines beginning with "#", or with all whitespace characters, are ignored
// and that the newline character and whitespace at the end is trimmed
// returns null on end of file, next line if one exists, and exits with error otherwise
char *getNextLine(int f) {
	char linetmp[9] = { 0 };

	if (f < 0 || f > MAX_OPEN_FILES-1) printErrorAndExit(4, "Invalid File Handle!\n");
	if (!curFile[f]) printErrorAndExit(4, "Cannot read next line: no file open!\n");

getnextline_goto:
	if (!fgets(curFileLine[f], MAX_LINE_LEN, curFile[f])) { closeFile(f); return NULL; }
	if (strlen(curFileLine[f]) < 1) goto getnextline_goto; // if empty, goto next line
	if (sscanf(curFileLine[f], " %8s", linetmp) < 1) goto getnextline_goto; // if all whitespace, goto next line
	if (linetmp[0] == '#') goto getnextline_goto; // if entire line is a perl comment, goto next line

	while (strlen(curFileLine[f]) > 0 && isspace(curFileLine[f][strlen(curFileLine[f])-1]))
		curFileLine[f][strlen(curFileLine[f])-1] = 0;

	if (strlen(curFileLine[f]) == 0) goto getnextline_goto;

	return curFileLine[f];
}

// END MINI FILE PARSING MODULE

