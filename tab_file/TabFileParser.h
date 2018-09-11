/*
 * TabFileParser.h
 *
 *  Created on: Mar 7, 2012
 *      Author: xinfeng
 */

#ifndef TABFILEPARSER_H_
#define TABFILEPARSER_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <stdint.h>
#include "TabGene.h"

namespace tab_file {

//template<typename Result>
class TabFileParser {
public:

	TabFileParser(const char* file, std::map<std::string,std::vector<TabGene> >& result);
	TabFileParser(std::string& file, std::map<std::string,std::vector<TabGene> >& result);
	virtual ~TabFileParser();
	void parse();

protected:
	void parseTokens(std::vector<std::string>& tokens);
	void skipHeadLines(size_t n, const char* prefix, std::istream& is);

private:
	TabFileParser();
	typedef size_t index_t;
	std::map<std::string,std::vector<TabGene> >& _result;
	std::string _file;
	static const index_t chr_col;
	static const index_t dir_col;
	static const index_t gene_start_col;
	static const index_t gene_end_col;
	static const index_t utr5_end_col;
	static const index_t utr3_start_col;
	static const index_t exon_start_col;
	static const index_t exon_end_col;
	static const index_t gene_name2;
};

} /* namespace tab_file */
#endif /* TABFILEPARSER_H_ */
