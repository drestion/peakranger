/*
 * sam_read.h
 *
 * based on BamAlignment from bamtools
 *  Created on: Jun 29, 2011
 *      Author: xin
 */

#ifndef SAM_READ_H_
#define SAM_READ_H_

#include <string>
#include <stdint.h>
#include "utils/exceptions.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>


class sam_read{

    // constructors & destructor
public:
    sam_read(void);
    void parseLine(std::string& line);
    ~sam_read(void);

    // Queries against alignment flags
public:
    bool IsDuplicate(void) const; // Returns true if this read is a PCR duplicate
    bool IsFailedQC(void) const; // Returns true if this read failed quality control
    bool IsFirstMate(void) const; // Returns true if alignment is first mate on read
    bool IsMapped(void) const; // Returns true if alignment is mapped
    bool IsMateMapped(void) const; // Returns true if alignment's mate is mapped
    bool IsMateReverseStrand(void) const; // Returns true if alignment's mate mapped to reverse strand
    bool IsPaired(void) const; // Returns true if alignment part of paired-end read
    bool IsPrimaryAlignment(void) const; // Returns true if reported position is primary alignment
    bool IsProperPair(void) const; // Returns true if alignment is part of read that satisfied paired-end resolution
    bool IsReverseStrand(void) const; // Returns true if alignment mapped to reverse strand
    bool IsSecondMate(void) const; // Returns true if alignment is second mate on read

    // Manipulate alignment flags
public:
    void SetIsDuplicate(bool ok); // Sets "PCR duplicate" flag
    void SetIsFailedQC(bool ok); // Sets "failed quality control" flag
    void SetIsFirstMate(bool ok); // Sets "alignment is first mate" flag
    void SetIsMateUnmapped(bool ok); // Sets "alignment's mate is mapped" flag
    void SetIsMateReverseStrand(bool ok); // Sets "alignment's mate mapped to reverse strand" flag
    void SetIsPaired(bool ok); // Sets "alignment part of paired-end read" flag
    void SetIsProperPair(bool ok); // Sets "alignment is part of read that satisfied paired-end resolution" flag
    void SetIsReverseStrand(bool ok); // Sets "alignment mapped to reverse strand" flag
    void SetIsSecondaryAlignment(bool ok); // Sets "position is primary alignment" flag
    void SetIsSecondMate(bool ok); // Sets "alignment is second mate on read" flag
    void SetIsUnmapped(bool ok); // Sets "alignment is mapped" flag

    // Tag data access methods
public:
    bool GetEditDistance(uint8_t& editDistance) const; // get "NM" tag data - contributed by Aaron Quinlan
    bool GetReadGroup(std::string& readGroup) const; // get "RG" tag data

    bool GetTag(const std::string& tag,
                std::string& destination);
    template<typename T> bool GetTag(const std::string& tag,
                                     T& destination);

    // Additional data access methods
public:
    int GetEndPosition(bool usePadded = false) const; // calculates alignment end position, based on starting position and CIGAR operations

    // 'internal' utility methods
private:
    static void SkipToNextTag(const char storageType,
                              char* &pTagData,
                              unsigned int& numBytesParsed);

    // Data members
public:
    std::string Name; // Read name
    int32_t Length; // Query length
    std::string QueryBases; // 'Original' sequence (as reported from sequencing machine)
    //    std::string AlignedBases; // 'Aligned' sequence (includes any indels, padding, clipping)
    //    std::string Qualities; // FASTQ qualities (ASCII characters, not numeric values)
    //    std::string TagData; // Tag data (accessor methods will pull the requested information out)
    std::string RefName; // ID number for reference sequence
    int32_t Position; // Position (0-based) where alignment starts
    //    uint16_t Bin; // Bin in BAM file where this alignment resides
    //    uint16_t MapQuality; // Mapping quality score
    uint32_t AlignmentFlag; // Alignment bit-flag - see Is<something>() methods to query this value, SetIs<something>() methods to manipulate
    //    std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
    int32_t MateRefID; // ID number for reference sequence where alignment's mate was aligned
    int32_t MatePosition; // Position (0-based) where alignment's mate starts
    //    int32_t InsertSize; // Mate-pair insert size


private:
    enum {
        sPAIRED = 1,
        sPROPER_PAIR = 2,
        sUNMAPPED = 4,
        sMATE_UNMAPPED = 8,
        sREVERSE = 16,

        sMATE_REVERSE = 32,
        sREAD_1 = 64,
        sREAD_2 = 128,
        sSECONDARY = 256,
        sQC_FAILED = 512,
        sDUPLICATE = 1024
    };
};
inline sam_read::sam_read(void) {
}


inline sam_read::~sam_read(void) {
}

// Queries against alignment flags
inline bool sam_read::IsDuplicate(void) const {
    return ((AlignmentFlag & sDUPLICATE) != 0);
}
inline bool sam_read::IsFailedQC(void) const {
    return ((AlignmentFlag & sQC_FAILED) != 0);
}
inline bool sam_read::IsFirstMate(void) const {
    return ((AlignmentFlag & sREAD_1) != 0);
}
inline bool sam_read::IsMapped(void) const {
    return ((AlignmentFlag & sUNMAPPED) == 0);
}
inline bool sam_read::IsMateMapped(void) const {
    return ((AlignmentFlag & sMATE_UNMAPPED) == 0);
}
inline bool sam_read::IsMateReverseStrand(void) const {
    return ((AlignmentFlag & sMATE_REVERSE) != 0);
}
inline bool sam_read::IsPaired(void) const {
    return ((AlignmentFlag & sPAIRED) != 0);
}
inline bool sam_read::IsPrimaryAlignment(void) const {
    return ((AlignmentFlag & sSECONDARY) == 0);
}
inline bool sam_read::IsProperPair(void) const {
    return ((AlignmentFlag & sPROPER_PAIR) != 0);
}
inline bool sam_read::IsReverseStrand(void) const {
    return ((AlignmentFlag & sREVERSE) != 0);
}
inline bool sam_read::IsSecondMate(void) const {
    return ((AlignmentFlag & sREAD_2) != 0);
}

// Manipulate alignment flags
inline void sam_read::SetIsDuplicate(bool ok) {
    if (ok) AlignmentFlag |= sDUPLICATE;
    else AlignmentFlag &= ~sDUPLICATE;
}
inline void sam_read::SetIsFailedQC(bool ok) {
    if (ok) AlignmentFlag |= sQC_FAILED;
    else AlignmentFlag &= ~sQC_FAILED;
}
inline void sam_read::SetIsFirstMate(bool ok) {
    if (ok) AlignmentFlag |= sREAD_1;
    else AlignmentFlag &= ~sREAD_1;
}
inline void sam_read::SetIsMateUnmapped(bool ok) {
    if (ok) AlignmentFlag |= sMATE_UNMAPPED;
    else AlignmentFlag &= ~sMATE_UNMAPPED;
}
inline void sam_read::SetIsMateReverseStrand(bool ok) {
    if (ok) AlignmentFlag |= sMATE_REVERSE;
    else AlignmentFlag &= ~sMATE_REVERSE;
}
inline void sam_read::SetIsPaired(bool ok) {
    if (ok) AlignmentFlag |= sPAIRED;
    else AlignmentFlag &= ~sPAIRED;
}
inline void sam_read::SetIsProperPair(bool ok) {
    if (ok) AlignmentFlag |= sPROPER_PAIR;
    else AlignmentFlag &= ~sPROPER_PAIR;
}
inline void sam_read::SetIsReverseStrand(bool ok) {
    if (ok) AlignmentFlag |= sREVERSE;
    else AlignmentFlag &= ~sREVERSE;
}
inline void sam_read::SetIsSecondaryAlignment(bool ok) {
    if (ok) AlignmentFlag |= sSECONDARY;
    else AlignmentFlag &= ~sSECONDARY;
}
inline void sam_read::SetIsSecondMate(bool ok) {
    if (ok) AlignmentFlag |= sREAD_2;
    else AlignmentFlag &= ~sREAD_2;
}
inline void sam_read::parseLine(std::string & readline) {
    uint32_t flag;
    size_t ploc;
    std::string chr, seq, line, chr2, seq2;
    uint32_t loc, loc2;
    line = readline;
    ploc = line.find('\t');
    if (ploc == std::string::npos) {
        throw DataLineNotValid(readline);
    }
    std::string name(line.begin(),
                line.begin() + ploc);
    std::string nameTrimmed(line.begin() + ploc,
                       line.end());
    std::stringstream iss(nameTrimmed);
    iss >> flag >> chr >> loc;
    iss >> seq;
    iss >> seq;
    iss >> seq;
    iss >> loc2;
    iss >> seq;
    iss >> seq;
    Name = name;
    Length = seq.size();
    QueryBases = seq;
    RefName = chr;
    Position = loc;
    AlignmentFlag = flag;
    MatePosition = loc2;
}

inline void sam_read::SetIsUnmapped(bool ok) {
    if (ok) AlignmentFlag |= sUNMAPPED;
    else AlignmentFlag &= ~sUNMAPPED;
}

#endif /* SAM_READ_H_ */
