/**
 * @file print_xml_doc.h
 *
 * @date Jun 13, 2011
 * @author ali
 */

#ifndef PRINT_XML_DOC_H_
#define PRINT_XML_DOC_H_

/**
 * This function reads xml files given in command line and prints them
 * to the stdout.
 *
 * @param argc number of command line args
 * @param argv the array of passed args
 */
void print_xml_doc(int argc, char* argv[]);
void dump_to_stdout(const char* pFilename);

#endif /* PRINT_XML_DOC_H_ */
