/*
 * cmd_option_parser.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: xin
 */

#include "cmd_option_parser.h"
using namespace std;
const string cmd_option_parser::format_bowtie = "bowtie";
const string cmd_option_parser::format_bed = "bed";
const string cmd_option_parser::format_sam = "sam";
const string cmd_option_parser::format_bam = "bam";
const string cmd_option_parser::format_eland = "eland";

string cmd_option_parser::getControl_dir() const {
    return _control_dir;
}

string cmd_option_parser::getControl_file() const {
    return _control_file;
}

double cmd_option_parser::getCut_off() const {
    return _p_cut_off;
}

double cmd_option_parser::getDelta() const {
    return _delta;
}

uint32_t cmd_option_parser::getExt_length() const {
    return _ext_length;
}

string cmd_option_parser::getFormat() const {
    return _format;
}

string cmd_option_parser::getMode() const {
    return _mode;
}

uint32_t cmd_option_parser::getNo_of_thread() const {
    return _no_of_thread;
}

bool cmd_option_parser::getNowig() const {
    return _nowig;
}

string cmd_option_parser::getOutput_dir() const {
    return _output_dir;
}

string cmd_option_parser::getOutput_file() const {
    return _output_file;
}

bool cmd_option_parser::getPad() const {
    return _pad;
}

string cmd_option_parser::getTreat_dir() const {
    return _treat_dir;
}

string cmd_option_parser::getTreat_file() const {
    return _treat_file;
}

uint32_t cmd_option_parser::getBandwidth() const {
    return _bandwidth;
}

bool cmd_option_parser::getHelp() const {
    return _helpRequested;
}

bool cmd_option_parser::getSpecified_output() const {
    return _outputSpecified;
}

bool cmd_option_parser::getUsing_chr_table() const {
    return _chrtableSpecified;
}

bool cmd_option_parser::getVersion() const {
    return _versionRequested;
}

void cmd_option_parser::setBandwidth(uint32_t _bandwidth) {
    this->_bandwidth = _bandwidth;
}

string cmd_option_parser::getChr_table_file() const {
    return _chr_table_file;
}

void cmd_option_parser::setChr_table_file(string _chr_table) {
    this->_chr_table_file = _chr_table;
    this->setChrtableSpecified(true);
}

void cmd_option_parser::setControl_dir(string _control_dir) {
    this->_control_dir = _control_dir;
}

void cmd_option_parser::setControl_file(string _control_file) {
    this->_control_file = _control_file;
}

void cmd_option_parser::setCut_off(double _p_cut_off) {
    this->_p_cut_off = _p_cut_off;
}

void cmd_option_parser::setDelta(double _delta) {
    this->_delta = _delta;
}

void cmd_option_parser::setExt_length(uint32_t _ext_length) {
    this->_ext_length = _ext_length;
}

void cmd_option_parser::setFormat(string _format) {
    this->_format = _format;
}

void cmd_option_parser::setHelp(bool _help) {
    this->_helpRequested = _help;
}

void cmd_option_parser::setMode(string _mode) {
    this->_mode = _mode;
}

void cmd_option_parser::setNo_of_thread(uint32_t _no_of_thread) {
    this->_no_of_thread = _no_of_thread;
}

void cmd_option_parser::setNowig(bool _nowig) {
    this->_nowig = _nowig;
}

void cmd_option_parser::setOutput_dir(string _output_dir) {
    this->_output_dir = _output_dir;
    this->setOutputSpecified(true);
}

void cmd_option_parser::setOutput_file(string _output_file) {
    this->_output_file = _output_file;
    this->setOutputSpecified(true);
}

void cmd_option_parser::setPad(bool _pad) {
    this->_pad = _pad;
}

void cmd_option_parser::setSpecified_output(bool _specified_output) {
    this->_outputSpecified = _specified_output;
}

void cmd_option_parser::setTreat_dir(string _treat_dir) {
    this->_treat_dir = _treat_dir;
}

void cmd_option_parser::setTreat_file(string _treat_file) {
    this->_treat_file = _treat_file;
}

void cmd_option_parser::setVersion(bool _version) {
    this->_versionRequested = _version;
}

bool cmd_option_parser::getChrtableSpecified() const {
    return _chrtableSpecified;
}

bool cmd_option_parser::getHelpRequested() const {
    return _helpRequested;
}

bool cmd_option_parser::getOutputSpecified() const {
    return _outputSpecified;
}

bool cmd_option_parser::getVersionRequested() const {
    return _versionRequested;
}

void cmd_option_parser::setChrtableSpecified(bool _chrtableSpecified) {
    this->_chrtableSpecified = _chrtableSpecified;
}

void cmd_option_parser::setHelpRequested(bool _helpRequested) {
    this->_helpRequested = _helpRequested;
}

void cmd_option_parser::setOutputSpecified(bool _outputSpecified) {
    this->_outputSpecified = _outputSpecified;
}

void cmd_option_parser::setVersionRequested(bool _versionRequested) {
    this->_versionRequested = _versionRequested;
}

bool cmd_option_parser::getVerboseRequested() const {
    return _verboseRequested;
}

void cmd_option_parser::setUsing_chr_table(bool _using_chr_table) {
    this->_chrtableSpecified = _using_chr_table;
}

void cmd_option_parser::setVerboseRequested(bool _verboseRequested) {
    this->_verboseRequested = _verboseRequested;
}

string cmd_option_parser::getTreatfilename() const {
    return _treatfilename;
}

void cmd_option_parser::setTreatfilename(string _treatfilename) {
    this->_treatfilename = _treatfilename;
}

string cmd_option_parser::getControlfilename() const {
    return _controlfilename;
}

void cmd_option_parser::setControlfilename(string _controlfilename) {
    this->_controlfilename = _controlfilename;
}

string cmd_option_parser::getControl_wig_file() const {
    return _control_wig_file;
}

string cmd_option_parser::getTreat_wig_file() const {
    return _treat_wig_file;
}

void cmd_option_parser::setControl_wig_file(string _control_wig_file) {
    this->_control_wig_file = _control_wig_file;
}

void cmd_option_parser::setTreat_wig_file(string _treat_wig_file) {
    this->_treat_wig_file = _treat_wig_file;
}

vector<string> cmd_option_parser::getChrs_to_parse() const {
    return _chrs_to_parse;
}

void cmd_option_parser::setChrs_to_parse(vector<string> _chrs_to_parse) {
    this->_chrs_to_parse = _chrs_to_parse;
}

void cmd_option_parser::setBinlength(uint32_t _binlength) {
    this->_binlength = _binlength;
}

string cmd_option_parser::getConfigFile() const {
    return _config_file;
}

void cmd_option_parser::setConfigFile(string _config_file) {
    this->_config_file = _config_file;
}

double cmd_option_parser::getFdrCutOff() const {
    return _fdr_cut_off;
}

void cmd_option_parser::setFdrCutOff(double _fdr_cut_off) {
    this->_fdr_cut_off = _fdr_cut_off;
}

string cmd_option_parser::getReportName() const {
    return _report_name;
}

string cmd_option_parser::getChr_table() const {
    return _chr_table_file;
}

void cmd_option_parser::setAc(int _ac) {
}

void cmd_option_parser::setAv(char **_av) {
}

void cmd_option_parser::setChr_table(string _chr_table) {
}

uint32_t cmd_option_parser::getBinlength() const {
    return _binlength;
}

void cmd_option_parser::setReportName(string _report_name) {
    this->_report_name = _report_name;
}

uint32_t cmd_option_parser::getHtmlRegionLength() const {
    return _html_region_length;
}

