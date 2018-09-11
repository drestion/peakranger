/*
 * chipseqhtmlreporter.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: xfeng
 */

#include "chipseqhtmlreporter.h"
#include "wiggle/wig.h"
#include "region_profile/profilezoom.h"
#include "tab_file/TabFileParser.h"
#include "tab_file/TabGene.h"
#include "tab_file/Region.h"
#include "utils/timer.h"
#include "utils/stl_helper.h"
#include "ggplay/HtmlAux.h"
#include "ggplay/ImgScriptRunner.h"
#include "concepts/RegionInt32.h"
#include <sys/stat.h> //todo: win
#include <sys/types.h>//todo: win
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <boost/algorithm/string.hpp>
using namespace utils;
using namespace std;
using namespace ggplay;
using namespace ranger::concepts;
using namespace tab_file;
chipseq_html_reporter::chipseq_html_reporter() {

    _binlength = 1000000; //for wig builder
    _d = 5000;
    _report_name = "reports";
    format_name();
}

void chipseq_html_reporter::print_head(ostream& os, uint32_t reportLength) {

    os
            << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">";
    os
            << "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=ISO-8859-1\">";

    os << "<title>Discovered peaks</title>";

    os << "<style type=\"text/css\">";
    os << "body {";
    os << "font-family: Arial, Helvetica, sans-serif;";
    os << "font-size: 14px;";
    os << "}";

    os << ".bold {";
    os << "font-weight: bold !important;";
    os << "}";

    os << ".pager {";
    os << "margin-left: auto;";
    os << "margin-right: auto;";
    os << "padding: 0px 100px 10px 20px;";
    os << "}";

    os << "/* define table skin */";
    os << "table.grid {";
    os << "margin: 0;";
    os << "padding: 0;";
    os << "border-collapse: separate;";
    os << "border-spacing: 0;";
    os << "width: 100%;";
    os << "}";

    os << "table.grid * {";
    os << "font: 11px Arial, Helvetica, sans-serif;";
    os << "vertical-align: top;";
    os << "text-align: left;";
    os << "}";

    os << "table.grid thead, table.grid .collapsible {";
    os << "background-color: #e6edc1;";
    os << "}";

    os << "table.grid th {";
    os << "color: #565770;";
    os << "padding: 4px 16px 4px 0;";
    os << "cursor: pointer;";
    os << "} ";

    os << "table.grid td {";
    os << "color: #565770;";
    os << "padding: 4px 6px;";
    os << "}";

    os << "table.grid th.headerSortUp {";
    os
            << "background: #e6edc1 url(./tablesorter_bkgrd.png) no-repeat 100% -80px;";
    os << "}";

    os << "table.grid th.headerSortUp span {";
    os << "background: #e6edc1 url(./tablesorter_bkgrd.png) no-repeat 0 -80px;";
    os << "}";

    os << "table.grid th.headerSortDown {";
    os
            << "background: #e6edc1 url(./tablesorter_bkgrd.png) no-repeat 100% 2px;";
    os << "}";

    os << "table.grid th.headerSortDown span {";
    os << "background: #e6edc1 url(./tablesorter_bkgrd.png) no-repeat 0 2px;";
    os << "}";

    os << "table.grid th span {";
    os << "padding: 4px 0 4px 6px;";
    os << "font-weight: bold;";
    os << "}";

    os << "table.grid a:link,";
    os << "table.grid a:visited,";
    os << "table.grid a:focus,";
    os << "table.grid a:hover {";
    os << "color: #565770;";
    os << "font-weight: bold;";
    os << "text-decoration: underline;";
    os << "}";

    os << "table.grid a:hover {";
    os << "color: #d07c57;";
    os << "}";

    os << "table.grid tr.even {";
    os << "background-color: #f0f0f0;";
    os << "}";

    os << "/* expand/collapse */";
    os << "table.grid .collapsible {";
    os << "padding: 0 0 3px 0;";
    os << "}";

    os << ".collapsible a.collapsed {";
    os << "margin: 2px;";
    os << "display: block;";
    os << "width: 15px;";
    os << "height: 15px;";
    os
            << "background: url(./scripts/tablesorter_expand.png) no-repeat 3px 3px;";
    os << "outline: 0;";
    os << "}";

    os << ".collapsible a.expanded {";
    os << "margin: 2px;";
    os << "display: block;";
    os << "width: 15px;";
    os << "height: 15px;";
    os
            << "background: url(./scripts/tablesorter_expand.png) no-repeat -24px 3px;";
    os << "outline: 0;";
    os << "}";
    os << "</style>";
    os
            << "<script type=\"text/javascript\" src=\"./scripts/jquery-1.2.3.js\"></script>";
    os
            << "<script type=\"text/javascript\" src=\"./scripts/jquery.tablesorter.mod.js\"></script>";
    os
            << "<script type=\"text/javascript\" src=\"./scripts/jquery.tablesorter.pager.js\"></script>";
    os
            << "<script type=\"text/javascript\" src=\"./scripts/jquery.tablesorter.collapsible.js\"></script>";
    os << "<script type=\"text/javascript\">";
    os << "$(document).ready(";
    os << "function (){";
    os << "$(\".tablesorter\")";
    os << ".collapsible(\"td.collapsible\", {";
    os << "collapse: true";
    os << "})";
    os << ".tablesorter({";
    os << "sortList: [[1,0],[2,0]],";
    os
            << "headers: {0: {sorter: false},1:{sorter:'text'},5: {sorter: 'scinot'}}";
    os << ", widgets: ['zebra']";
    os << ", onRenderHeader: function (){";
    os << "this.wrapInner(\"<span></span>\");";
    os << "}";
    os << ", debug: false";
    os << "})";
    os << ".tablesorterPager({container: $(\"#pager\"), positionFixed: false})";
    os << ";";
    os << "}";
    os << ");";
    os << "</script>";
    os << "</head>";
    os << "<body>";
    os << "<h1>";
    os << "Discovered peaks";
    os << "</h1>";
    printCitation(os);
    os << "<div id=\"pager\" class=\"pager\" >";
    os << "<form>";
    os << "<input type=\"button\" value=\"&lt;&lt;\" class=\"first\">";
    os << "<input type=\"button\" value=\"&lt;\" class=\"prev\">";
    os << "<input type=\"text\" class=\"pagedisplay\">";
    os << "<input type=\"button\" value=\"&gt;\" class=\"next\">";
    os << "<input type=\"button\" value=\"&gt;&gt;\" class=\"last\">";
    os << "<select class=\"pagesize\">";
    os << "<option selected=\"selected\" value=\"10\">10</option>";
    os << "<option value=\"50\">50</option>";
    os << "<option value=\"100\">100</option>";
    os << "<!--<option value=\"40\">40</option> </-->";
    os << "</select>";
    os << "</form>";
    os << "</div>";

    os
            << "<table class=\"grid tablesorter\" cellspacing=\"0\" style=\"width: auto;\">";
    os << "<colgroup>";
    os << "<col width=\"19\">";
    os << "<col width=\"30\">"; //chr
    os << "<col width=\"90\">"; //start
    os << "<col width=\"90\">"; //end
    os << "<col width=\"180\">"; //gene name 3k
//    os << "<col width=\"180\">"; //gene name 10k
    os << "<col width=\"180\">"; //summits
    os << "<col width=\"90\">"; //p
    os << "<col width=\"90\">"; //q
    os << "<col width=\"45\">"; //reads
    os << "<col width=\"45\">"; //control reads
    os << "</colgroup>";

    os << "<thead>";
    os << "<tr>";
    os << "<th style=\"padding: 0; margin: 0;\" class=\"\"></th>";
    os << "<th class=\"header\"><span>";
    os << "Chr";
    os << "</span></th>";
    os << "<th class=\"header\"><span>";
    os << "Start";
    os << "</span></th>";
    os << "<th class=\"header\"><span>";
    os << "End";
    os << "</span></th>";

    os << "<th class=\"header\"><span>";
    os << "Genes(" << "within " << reportLength / 1000 << "kbp)";
    os << "</span></th>";

    os << "<th class=\"header\"><span>";
    os << "Summits";
    os << "</span></th>";
    os << "<th class=\"header headerSortDown\"><span>";
    os << "P-value";
    os << "</span></th>";
    os << "<th class=\"header\"><span>";
    os << "FDR";
    os << "</span></th>";
    os << "<th class=\"header\"><span>";
    os << "Reads1";
    os << "</span></th>";
    os << "<th class=\"header\"><span>";
    os << "Reads2";
    os << "</span></th>";
    os << "</tr>";
    os << "</thead>";
    os << "<tbody>";
}

void chipseq_html_reporter::print_r_head(ostream& os) {

    os << "drawGene <- "
            << "function(xl,y,y2,utr5,utr3,e_s,e_e,gene_id,strand=TRUE,trs_h=4){\n";
    os << "xl1=xl[1];\n";
    os << "xl2=xl[2];\n";
    os << "e_col=rgb(51/255,51/255,51/255)\n";
    os << "\n";
    os << "#intron 5'\n";
    os << "if_s=e_e[1:length(e_e)-1]+1;\n";
    os << "if_e=e_s[2:length(e_s)]-1;\n";
    os << "\n";
    os << "\n";
    os << "##plot gene model\n";
    os << "# trs_h = 4\n";
    os << "# par(fig=c(0.05,0.9,y,y2),new=TRUE,mgp=c(0.1,0.1,0.1))\n";
    os << "   ###plot utr5\n";
    os << "       lines(utr5,\n";
    os << "           rep(trs_h,2),\n";
    os << "           xlim=c(xl1,xl2),\n";
    os << "           lwd=6,\n";
    os << "           lend=2,\n";
    os << "           ylim=c(0,10),\n";
    os << "           bty=\"n\",\n";
    os << "           type=\"l\",\n";
    os << "           yaxt=\"n\" ,xaxt=\"n\",\n";
    os << "           col=e_col,\n";
    os << "           xlab=\"\",\n";
    os << "           ylab=\"\")\n";
    os << "\n";
    os << "   ###plot introns\n";
    os << "   for(i in 1:length(if_s)){\n";
    os << "        lines(c(if_s[i],if_e[i]),\n";
    os << "       rep(trs_h,2),\n";
    os << "       xlim=c(xl1,xl2),\n";
    os << "       lwd=1,\n";
    os << "       lend=2,\n";
    os << "       ylim=c(0,10),\n";
    os << "       bty=\"n\",\n";
    os << "       type=\"l\",\n";
    os << "       yaxt=\"n\" ,xaxt=\"n\",\n";
    os << "       col=e_col,\n";
    os << "       xlab=\"\",\n";
    os << "       ylab=\"\")\n";
    os << "       }\n";
    os << "       \n";
    os << "   ###plot exons\n";
    os << "   for(i in 1:length(e_s)){\n";
    os << "       lines(c(e_s[i],e_e[i]),\n";
    os << "           rep(trs_h,2),\n";
    os << "           xlim=c(xl1,xl2),\n";
    os << "           lwd=10,\n";
    os << "           lend=2,\n";
    os << "           ylim=c(0,10),\n";
    os << "           bty=\"n\",\n";
    os << "           type=\"l\",\n";
    os << "           yaxt=\"n\" ,xaxt=\"n\",\n";
    os << "           col=e_col,\n";
    os << "           xlab=\"\",\n";
    os << "           ylab=\"\")\n";
    os << "       }\n";
    os << "       \n";
    os << "   ###plot utr3\n";
    os << "       lines(utr3,\n";
    os << "           rep(trs_h,2),\n";
    os << "           xlim=c(xl1,xl2),\n";
    os << "           lwd=6,\n";
    os << "           lend=2,\n";
    os << "           ylim=c(0,10),\n";
    os << "           bty=\"n\",\n";
    os << "           type=\"l\",\n";
    os << "           yaxt=\"n\" ,xaxt=\"n\",\n";
    os << "           col=e_col,\n";
    os << "           xlab=\"\",\n";
    os << "           ylab=\"\")\n";
    os << "   ##plot gene orientation \n";
    os << "   ar_ratio=0.015;\n";
    os
            << "   ar_xs=seq(max(if_s[1],e_s[1]),max(if_e[length(if_e)],e_e[length(e_e)]),by=ar_ratio*(xl2-xl1))\n";

    os << "   ar_xe=ar_xs+20\n";
    os << "   if(strand==FALSE){\n";
    os << "   ar_xe=ar_xs-20;}\n";
    os << "   ar_ys=rep(trs_h,length(ar_xs))\n";
    os << "   ar_ye=ar_ys\n";
    os << "   \n";
    os << "   arrows(ar_xs,\n";
    os << "       ar_ys,\n";
    os << "       ar_xe,\n";
    os << "       ar_ye,\n";
    os << "       length=0.05,\n";
    os << "       col=\"grey31\",\n";
    os << "   )\n";
    os << "   \n";
    os << "   ##plot gene name\n";
    os << "   \n";
    os
            << "txl=max(c(utr5[1],xl1));\ntxr=min(c(utr3[2],xl2));\nloc_text=mean(c(txl,txr));\n";
//    os << "   if(xl1 >= utr5[1] && xl2 <= utr3[2] ){\n";
//    os << "       loc_text =mean(xl);\n";
//    os << "   }\n";
//    os << "   else if(xl1 > utr5[1] && utr3[2] < xl2){\n";
//    os << "       loc_text=utr3[2];\n";
//    os << "   }\n";
//    os << "   else if(xl1 < utr5[1] && utr3[2] < xl2){\n";
//    os << "       loc_text =mean(e_s);\n";
//    os << "   }\n";
//    os << "   else if(xl1 < utr5[1] && xl2 < utr3[2]){\n";
//    os << "       loc_text=utr5[1];\n";
//    os << "   }\n";
    os << "   text(loc_text,1,gene_id,cex=0.6);\n";
    os << "}\n";
    os << "\n";
    os << " wid=800;";
    os << " hei=wid*0.618;";
}

void chipseq_html_reporter::print_r_tc_drawing(ostream& os) {
    os << "width=wid,height=hei,\n";
    os << "units=\"px\",\n";
    os << "pointsize=25);\n";
    os << "\n";
    os << "par(mar=c(0.5,0,0,0),oma=rep(0,4));\n";
    os << "plot.new();\n";
    os << "xl=c(mean(pk)-d,mean(pk)+d);\n";
    os << "ga=c(bb,cc);\n";
    os << "names<-c(\"Treatment\",\"Control\");\n";
    os << "colors<-c(\n";
    os << "\"#CD0000\",\n";
    os << "\"#7d7d7d\",\n";
    os << "\"#6c7b8b\",\n";
    os << "\"#bc8f8f\",\n";
    os << "\"#C51B8A\",\n";
    os << "\"#7A0177\");\n";
    os << "col_ax=\"dimgray\";\n";
    os << "n<-0;\n";
    os << "for(i in seq(1,(length(ga)),2)){\n";
    os << "   n<-n+1;\n";
    os << "   y<-0.95-n*0.4;\n";
    os << "   y2<-1-0.25*n;\n";
    os << "   if(i < 2){\n";
    os << "       y = 0.65-0.1;\n";
    os << "       y2 = 0.85-0.1;\n";
    os << "   }\n";
    os << "   else{\n";
    os << "       y = 0.25+0.08;\n";
    os << "       y2 = 0.45+0.08;\n";
    os << "   }\n";
    os << "   par(fig=c(0.05,0.9,y,y2), \n";
    os << "   mar=c(0.22,0,0,0),\n";
    os << "   mgp=c(0.1,0.1,0.1),\n";
    os << "   new=TRUE);\n";
    os << "   xx=unlist(ga[i],use.names=F);\n";
    os << "   yy=unlist(ga[i+1],use.names=F);\n";
    os << "   xx_ind=(xx>=xl[1] & xx<=xl[2]);\n";
    os << "   if( i < 2){\n";
    os << "       yl=c(0,max(yy[xx_ind]));\n";
    os << "   }\n";
    os
            << "   plot(xx[xx_ind],yy[xx_ind],bty=\"n\",type=\"h\",xlim=xl,ylim=yl,col=colors[n],\n";
    os << "   yaxt=\"n\" ,\n";
    os << "   xaxt=\"n\",\n";
    os << "   xlab=\"\",\n";
    os << "   ylab=\"\");\n";
    os << "   ti.x<-yl;\n";
    os << "   ti.l<-ti.x;\n";
    os
            << "   axis(2, at=ti.x,labels=ti.l,tck=-0.1, col.axis=col_ax, las=0,cex.axis=0.65);\n";
    os << "   mtext(names[n],EAST<-4,las=2,padj=1,cex=0.55,col=\"gray17\");\n";
    os << "}\n";
    os << "##plot the peak and the gene annotation\n";
    os << "y = y-0.245;\n";
    os << "y2 = y2-0.245;\n";
    os
            << "par(fig=c(0.05,0.9,y,y2),mgp=c(0.1,0.1,0.1), mar=c(0.22,0,0,0),new=TRUE);\n";
    os << "\n";
    os << "xl1=xl[1];\n";
    os << "xl2=xl[2];\n";
    os << "trs_h = 4;\n";
    os << "plot(c(xl1,xl2),\n";
    os << "       rep(trs_h,2),\n";
    os << "       xlim=c(xl1,xl2),\n";
    os << "       lwd=0,\n";
    os << "       lend=2,\n";
    os << "       ylim=c(0,10),\n";
    os << "       bty=\"n\",\n";
    os << "       type=\"p\",\n";
    os << "       yaxt=\"n\" ,xaxt=\"n\",\n";
    os << "       col=\"white\",\n";
    os << "       xlab=\"\",\n";
    os << "       ylab=\"\");\n";
}

void chipseq_html_reporter::filter_in_range_anno(
        std::vector<tab_file::TabGene> & genes_to_filter,
        std::vector<tab_file::TabGene> & genes_to_print) {
    sort(genes_to_filter.begin(), genes_to_filter.end(), tab_file::name_lex_lt);
    unique_copy(genes_to_filter.begin(), genes_to_filter.end(),
            std::back_inserter(genes_to_print), tab_file::sameID);
}

void chipseq_html_reporter::addCitation(const std::string& citation) {
    _cite = citation;
}

void chipseq_html_reporter::printCitation(std::ostream& os) {
    os << "<p>For further information on PeakRanger please visit ";
    os << "<a href=\"http://ranger.sourceforge.net/\"";
    os << ">http://ranger.sourceforge.net</a>.</p>";
    os << "<p>If you use PeakRanger in your research, please cite:<br>";
    os << _cite << "      </p>";
}

void chipseq_html_reporter::get_genes_to_filter(const string & chr,
        tab_file::Region & pkrg, vector<tab_file::TabGene> & genes_to_filter) {
    pair<vector<TabGene>::iterator, vector<TabGene>::iterator> bounds;
    TabGene tgg;
    tgg.utr5.setL(pkrg.getL());
    tgg.utr3.setR(pkrg.getR());
    bounds = equal_range(_genes[chr].begin(), _genes[chr].end(), tgg);
    copy(bounds.first, bounds.second, std::back_inserter(genes_to_filter));
}

void chipseq_html_reporter::print_gene_annot(const called_peak& pk,
        const string& chr, ostream& os) {
    size_t mid = (pk.first + pk.second) / 2;
    size_t l = 0;
    if (mid > _d) {
        l = mid - _d;
    }
    tab_file::Region pkrg(l, mid + _d);
    size_t i;
    vector<TabGene> genes_to_filter, genes_to_print;
    get_genes_to_filter(chr, pkrg, genes_to_filter);
    filter_in_range_anno(genes_to_filter, genes_to_print);

    vector<TabGene>::iterator itt = genes_to_print.begin();
    for (; itt != genes_to_print.end(); itt++) {
        os.flush();
        os << "e_id=\"";
        os << itt->name;
        os << "\";\n";
        os << "e_s=c(";
        i = 0;
        while (i < itt->exons.size() - 1) {
            os << itt->exons[i++].getL() << ", ";
        }
        os << itt->exons[i].getL();
        os << ");\n";

        os << "e_e=c(";
        i = 0;
        while (i < itt->exons.size() - 1) {
            os << itt->exons[i++].getR() << ", ";
        }
        os << itt->exons[i].getR();
        os << ");\n";

        os << "utr3=c(";
        os << itt->utr3;
        os << ");\n";
        os << "utr5=c(";
        os << itt->utr5;
        os << ");\n";
        if (itt->dir) {
            os << "drawGene(xl,y,y2,utr5,utr3,e_s,e_e,e_id,TRUE);\n";
        } else {
            os << "drawGene(xl,y,y2,utr5,utr3,e_s,e_e,e_id,FALSE);\n";
        }
    }
}
void chipseq_html_reporter::print_axis_to_end(ostream& os) {
    os
            << "lines(pk, rep(9,2),lwd=10,col=rgb(30/255,144/255,255/255),lend=1);\n";
    os << "cord=xl[1]\n";
    os
            << "ti.x<-c(xl[1],mean(c(mean(xl),xl[1])), mean(xl),mean(c(mean(xl),xl[2])),xl[2])\n";
    os << "ti.l<-ti.x;\n";
    os
            << "axis(1, at=ti.x,labels=ti.l, col.axis=col_ax, las=0,cex.axis=0.5,tck=-0.1)\n";
    os << "dev.off();\n";
}

void chipseq_html_reporter::print_css(cmd_option_parser& option) {
    string dir = option.getOutput_dir();
    string dir_r = dir + "/" + _report_name;
    string scripts = dir_r + "/scripts/";
    string jq = scripts + "jquery-1.2.3.js";
    string md = scripts + "jquery.tablesorter.mod.js";
    string pg = scripts + "jquery.tablesorter.pager.js";
    string cl = scripts + "jquery.tablesorter.collapsible.js";

    ofstream jqo(jq.c_str()), mdo(md.c_str()), pgo(pg.c_str()), clo(cl.c_str());

    ggplay::print_css_col(clo);
    ggplay::printJQTabmod(mdo);
    ggplay::printJQTabmodpager(pgo);
    ggplay::printJQuery(jqo);

}
void chipseq_html_reporter::print_peak(called_peak & pk, string chr,
        ostream& os) {

    std::ostringstream pkname;
    std::stringstream geneNamess;
    pkname << "./imgs/" << chr << "_" << pk.first << "_" << pk.second;
    os << " <tr class=\"odd\"> "
            "<td rowspan=\"2\"  class=\"collapsible\"> "
            "<a href=\"\" class=\"collapsed\"></a></td>"
            " <td rowspan=\"2\"  class=\"collapsible_alt\">";

    os << chr;
    os << "</td><td  >";
    os << pk.first;
    os << "</td><td  >";
    os << pk.second;
    os << "</td><td  >";
    RegionUint32 rg(pk.first, pk.second);
    vector<TabGene> overlapped;
    _nbgf.getOverlappedGenes(chr, rg, overlapped);
    foreach(TabGene& g, overlapped) {
        geneNamess << g.name << ",";
    }
    string geneNames(geneNamess.str());
    boost::replace_last(geneNames, ",", "");
    os << geneNames;
    os << "</td><td>";
    if (pk.summits.size() > 1) {
        std::copy(pk.summits.begin(), pk.summits.end(),
                ostream_iterator<uint32_t>(os, ","));
    } else if (pk.summits.size() > 0) {
        os << pk.summits.at(0);
    }
    os << "</td><td>";
    os << pk.p;
    os << "</td><td>";
    os << pk.q;
    os << "</td><td>";
    os << pk.treads;
    os << "</td><td align=\"center\">";
    os << pk.creads;

    os << "</td></tr><tr class=\"expand-child odd\"> "
            << "<td colspan=\"4\" style=\"display: none; \"> "
            << "<div class=\"bold\">Details of the peak</div><div><img src=\"";
    os << pkname.str() << ".png\" alt=\"" << pkname.str() << "\" /> </div>"
            "<a href=\"" << pkname.str() << ".png\" target=\"_blank\">"
            << "<p style=\"font-weight: bold;\">Download the image</p></a></td></tr>\n";
}

void chipseq_html_reporter::prepare_wigs(Reads& treads, Reads& creads,
        string& chr, cmd_option_parser& option, wigs& wt, wigs& wc) {
    _wig._binned_wig_compiler(_binlength, treads.getReadlength(),
            option.getExt_length(), treads.pos_reads.begin_of(chr),
            treads.pos_reads.end_of(chr), treads.neg_reads.begin_of(chr),
            treads.neg_reads.end_of(chr), wt);

    _wig._binned_wig_compiler(_binlength, creads.getReadlength(),
            option.getExt_length(), creads.pos_reads.begin_of(chr),
            creads.pos_reads.end_of(chr), creads.neg_reads.begin_of(chr),
            creads.neg_reads.end_of(chr), wc);
}

void chipseq_html_reporter::print_end(ostream& os) {
    os << "</tbody>";
    os << "</table>";
    os << "</body></html>";
}

void chipseq_html_reporter::prepare_dir(cmd_option_parser& option,
        const string& reportname) {
    //todo: win
    string dir = option.getOutput_dir();
    string dir_r = dir + "/" + reportname;
    mkdir(dir_r.c_str(), 0777);
    string scripts = dir_r + "/scripts";
    mkdir(scripts.c_str(), 0777);
    string imgs = dir_r + "/imgs";
    mkdir(imgs.c_str(), 0777);
}

void chipseq_html_reporter::run_pk_img_script(string& rf) {
    //todo: win32

    string R = "R CMD BATCH " + rf + " /dev/null";
    int r1 = system(R.c_str());
    string rm = "rm -f " + rf;
    r1 = system(rm.c_str());
}

void chipseq_html_reporter::run_img_scripts(cmd_option_parser& option,
        map<string, vector<called_peak> > & peaks) {
    string dir = option.getOutput_dir();
    string imgs = dir + "/" + _report_name + "/imgs/";
    map<string, vector<called_peak> >::iterator it;
    it = peaks.begin();
    if (option.getVerboseRequested()) {
        cout << "Generating details of peaks...\n";
        cout.flush();
    }
    for (; it != peaks.end(); it++) {
        string chr = it->first;
        vector<called_peak>::iterator it;
        it = peaks[chr].begin();
        vector<string> imgScripts;
        for (; it != peaks[chr].end(); it++) {
            std::ostringstream pkname;
            pkname << imgs << chr << "_" << it->first << "_" << it->second;
            string rs = pkname.str();
            imgScripts.push_back(rs);
        }

        ImgScriptRunner runner(imgScripts, option.getNo_of_thread());
        runner.run();
        if (option.getVerboseRequested()) {
            cout << "Finished HTML reports for "<<chr<<"\n";
        }
    }
    if (option.getVerboseRequested()) {
        cout << "HTML reports complete\n";
    }
}

void chipseq_html_reporter::print_html_table(
        map<string, vector<called_peak> > & peaks, cmd_option_parser & option,
        ostream& os) {
    map<string, vector<called_peak> >::iterator it;
    it = peaks.begin();
    for (; it != peaks.end(); it++) {
        string chr = it->first;
        print_chr(chr, peaks, os);
    }
}

void chipseq_html_reporter::print_chr(string & chr,
        map<string, vector<called_peak> > & peaks, ostream& os) {
    vector<called_peak>::iterator it;
    it = peaks[chr].begin();
    for (; it != peaks[chr].end(); it++) {
        print_peak(*it, chr, os);
    }
}

void chipseq_html_reporter::print_pk_img_wig(wigs& wt, wigs& wc,
        called_peak& pk, ostream& os) {
    /*  bb=data.frame(V1=c(),V2=c())
     cc=data.frame(V1=c(),V2=c())
     */
    ostringstream oss;
    os << "bb=data.frame(V1=c(";
    wigs::iterator it;
    uint32_t _xl1, _xl2, m;
    m = ((pk.first + pk.second) / 2);
    if (m > (_d)) {
        _xl1 = m - (_d);
    } else {
        _xl1 = 0;
    }
    _xl2 = m + (_d);
    it = wt.begin();
    for (; it != wt.end(); it++) {
        if (it->getP() >= _xl1 && it->getP() <= _xl2) {
            oss << it->getP() << ",";
        }
    }
    oss << "),V2=c(";
    it = wt.begin();
    for (; it != wt.end(); it++) {
        if (it->getP() >= _xl1 && it->getP() <= _xl2) {
            oss << it->getS() << ",";
        }
    }
    oss << "));";
    oss << "cc=data.frame(V1=c(";
    it = wc.begin();
    for (; it != wc.end(); it++) {
        if (it->getP() >= _xl1 && it->getP() <= _xl2) {
            oss << it->getP() << ",";
        }
    }
    oss << "),V2=c(";
    it = wc.begin();
    for (; it != wc.end(); it++) {
        if (it->getP() >= _xl1 && it->getP() <= _xl2) {
            oss << it->getS() << ",";
        }
    }
    oss << "));";
    string final = oss.str();
    boost::replace_all(final, ",),", "),");
    boost::replace_all(final, ",))", "))");
    os << final;
}
void chipseq_html_reporter::print_pk_img_para(ostream& os, called_peak& pk,
        string& rf) {
    os << "pk =c(" << pk.first << "," << pk.second << ");";
    os << "d = " << _d << ";";
    os << "png(\"" << rf;
    os << "\",";
}
void chipseq_html_reporter::print_peak_img_script(called_peak & pk, string rf,
        wigs& wt, wigs& wc, ostream& os, string chr) {
    print_r_head(os);
    print_pk_img_wig(wt, wc, pk, os);
    print_pk_img_para(os, pk, rf);
    print_r_tc_drawing(os);
    print_gene_annot(pk, chr, os);
    print_axis_to_end(os);
}

void chipseq_html_reporter::print_img_script(Reads& treads, Reads& creads,
        map<string, vector<called_peak> > & peaks, cmd_option_parser& option) {

    wigs wt, wc;
    string dir = option.getOutput_dir();
    string imgs = dir + "/" + _report_name + "/imgs/";

    map<string, vector<called_peak> >::iterator it;
    it = peaks.begin();
    for (; it != peaks.end(); it++) {
        string chr = it->first;
        wt.clear();
        wc.clear();
        prepare_wigs(treads, creads, chr, option, wt, wc);
        vector<called_peak>::iterator it;
        it = peaks[chr].begin();
        for (; it != peaks[chr].end(); it++) {
            std::ostringstream pkname;
            pkname << imgs << chr << "_" << it->first << "_" << it->second;
            ofstream os(pkname.str().c_str());
            pkname << ".png";
            print_peak_img_script(*it, pkname.str(), wt, wc, os, chr);
        }
    }
}

void chipseq_html_reporter::prepare_report_name(cmd_option_parser & option) {
    string date;
    getDate(date);
    _report_name = option.getTreatfilename() + "_report_" + date;
    boost::to_upper(_report_name);
    boost::replace_all(_report_name, ".", "_");
    boost::replace_all(_report_name, " ", "_");
}

void chipseq_html_reporter::generate_report(Reads & treads, Reads & creads,
        map<string, vector<called_peak> > & peaks, cmd_option_parser & option) {
    prepare_report_name(option);
    prepare_dir(option, _report_name);
    string dir = option.getOutput_dir();
    string dir_r = dir + "/" + _report_name;
    string idx = dir_r + "/index.html";
    string exp_img = dir_r + "/scripts/tablesorter_expand.png";
    _gene_anno_file = option.getGeneAnnoFile();
    setRegionLength(option.getHtmlRegionLength());
    ofstream os(idx.c_str());
    ofstream eimg(exp_img.c_str(), ios::binary);
    _nbgf.setAnnoFile(option.getGeneAnnoFile());
    _nbgf.setSearchSpan(option.getHtmlRegionLength());
    ggplay::print_exp_img(eimg);
    print_head(os, _d);
    print_html_table(peaks, option, os);
    print_end(os);
    print_css(option);
    prepare_gene_annot();
    print_img_script(treads, creads, peaks, option);
    run_img_scripts(option, peaks);
}

void chipseq_html_reporter::format_name() {
    boost::replace_all(_report_name, " ", "_");
    boost::replace_all(_report_name, "\t", "_");
}

void chipseq_html_reporter::setRegionLength(uint32_t _d) {
    this->_d = _d;
}

void chipseq_html_reporter::prepare_gene_annot() {
    tab_file::TabFileParser tabparser(_gene_anno_file, _genes);
    tabparser.parse();
}
