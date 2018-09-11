/*
 * cmd_option_parser.h
 *
 *  Created on: Jun 6, 2011
 *      Author: xin
 */

#ifndef CMD_OPTION_PARSER_H_
#define CMD_OPTION_PARSER_H_

#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
namespace po = boost::program_options;

namespace boost {
namespace program_options {
typedef options_description opt;
typedef positional_options_description p_opt;
}
}
//@to_be_depreciated
class cmd_option_parser {
public:
    cmd_option_parser() {

    }
    cmd_option_parser(int argc, char** argv) :
            _pad(false), _nowig(true), _verboseRequested(false), _helpRequested(
                    false), _versionRequested(false), _chrtableSpecified(false), _outputSpecified(
                    false), _estimateExt(false), _ac(argc), _av(argv) {
        maxThreads = boost::thread::hardware_concurrency();
    }

    virtual ~cmd_option_parser() {
    }

    virtual void print_option_file(std::ostream& os) const {
    }

    virtual void printHelp() const {

    }

    virtual void parse() = 0;
public:
    uint32_t fragmentSize;
    uint32_t slidingWinSize;
    uint32_t movingStep;
    uint32_t minCount;
    bool isStrandSensitiveMode;
    uint32_t outputNum;
    int randomSeed;
    double minScore;
    uint32_t bootstrapPass;
    double smoothingFactor;
public:
    std::string getChr_table() const;
    std::string getControl_dir() const;
    std::string getControl_file() const;
    double getCut_off() const;
    double getDelta() const;
    uint32_t getExt_length() const;
    std::string getFormat() const;
    std::string getMode() const;
    uint32_t getNo_of_thread() const;
    bool getNowig() const;
    std::string getOutput_dir() const;
    std::string getOutput_file() const;
    bool getPad() const;
    std::string getTreat_dir() const;
    std::string getTreat_file() const;

    uint32_t getBandwidth() const;

    bool getHelp() const;
    bool getSpecified_output() const;
    bool getUsing_chr_table() const;
    bool getVersion() const;
    void setAc(int _ac);
    void setAv(char **_av);
    void setBandwidth(uint32_t _bandwidth);
    void setChr_table(std::string _chr_table);
    void setControl_dir(std::string _control_dir);
    void setControl_file(std::string _control_file);
    void setCut_off(double _p_cut_off);
    void setDelta(double _delta);
    void setExt_length(uint32_t _ext_length);
    void setFormat(std::string _format);
    void setHelp(bool _help);
    void setMode(std::string _mode);
    void setNo_of_thread(uint32_t _no_of_thread);
    void setNowig(bool _nowig);
    void setOutput_dir(std::string _output_dir);
    void setOutput_file(std::string _output_file);
    void setPad(bool _pad);
    void setSpecified_output(bool _specified_output);
    void setTreat_dir(std::string _treat_dir);
    void setTreat_file(std::string _treat_file);
    void setUsing_chr_table(bool _using_chr_table);

    void setVersion(bool _version);
    bool getChrtableSpecified() const;
    bool getHelpRequested() const;
    bool getOutputSpecified() const;
    bool getVersionRequested() const;
    void setChrtableSpecified(bool _chrtableSpecified);
    void setHelpRequested(bool _helpRequested);
    void setOutputSpecified(bool _outputSpecified);
    void setVersionRequested(bool _versionRequested);
    std::string getChr_table_file() const;
    void setChr_table_file(std::string _chr_table_file);
    bool getVerboseRequested() const;
    void setVerboseRequested(bool _verboseRequested);
    std::string getTreatfilename() const;
    void setTreatfilename(std::string _treatfilename);
    std::string getControlfilename() const;
    void setControlfilename(std::string _controlfilename);
    std::string getControl_wig_file() const;
    std::string getTreat_wig_file() const;
    void setControl_wig_file(std::string _control_wig_file);
    void setTreat_wig_file(std::string _treat_wig_file);
    std::vector<std::string> getChrs_to_parse() const;
    void setChrs_to_parse(std::vector<std::string> _chrs_to_parse);
    uint32_t getBinlength() const;
    void setBinlength(uint32_t _binlength);
    std::string getConfigFile() const;
    void setConfigFile(std::string _config_file);
    double getFdrCutOff() const;
    void setFdrCutOff(double _fdr_cut_off);
    std::string getReportName() const;
    void setReportName(std::string _report_name);
    uint32_t getHtmlRegionLength() const;
    void setHtmlRegionLength(uint32_t _html_region_length);

    std::string getGeneAnnoFile() const {
        return _gene_anno_file;
    }

    void setGeneAnnoFile(std::string geneAnnoFile) {
        _gene_anno_file = geneAnnoFile;
    }
    bool getEstimateExt() const {
        return _estimateExt;
    }

    void setEstimateExt(bool _estimateExt) {
        this->_estimateExt = _estimateExt;
    }

    uint32_t getPeakHeightCutoff() const {
        return _threshold;
    }
    void setPeakHeightCutoff(uint32_t cut) {
        this->_threshold = cut;
    }
public:
    static const std::string format_bowtie;
    static const std::string format_bed;
    static const std::string format_sam;
    static const std::string format_bam;
    static const std::string format_eland;

protected:

    std::string _treat_file;
    std::string _treat_dir;
    std::string _treatfilename;
    std::string _control_file;
    std::string _control_dir;
    std::string _controlfilename;
    std::string _output_file;
    std::string _output_dir;
    std::string _control_wig_file;
    std::string _config_file;
    std::string _treat_wig_file;
    std::string _chr_table_file;
    std::string _format;
    std::string _mode;
    std::string _split_str;
    std::string _verbose_str;
    std::string _pad_str;
    std::string _nowig_str;
    std::string _report_name;
    std::string _gene_anno_file;

    double _p_cut_off;
    double _delta;
    double _fdr_cut_off;
    uint32_t _no_of_thread;
    uint32_t _ext_length;
    uint32_t _bandwidth;
    uint32_t _threshold;
    uint32_t _binlength;
    uint32_t _html_region_length;
    bool _pad;
    bool _nowig;
    bool _verboseRequested;
    bool _helpRequested;
    bool _versionRequested;
    bool _chrtableSpecified;
    bool _outputSpecified;
    bool _estimateExt;
    int _ac;
    char** _av;
    std::vector<std::string> _chrs_to_parse;
    uint32_t maxThreads;
};

#endif /* CMD_OPTION_PARSER_H_ */
