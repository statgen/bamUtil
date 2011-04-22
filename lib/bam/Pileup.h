/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __PILEUP_H__
#define __PILEUP_H__

template <class PILEUP_TYPE>
class defaultPileup
{
public:
    bool operator() (PILEUP_TYPE& element)
    {
        element.analyze();
        return(true);
    }
    void analyze(PILEUP_TYPE element)
    {
        element.analyze();
    }
};

template <class PILEUP_TYPE>
void defaultPileupAnalyze(PILEUP_TYPE& element)
{
    element.analyze();
}


template <class PILEUP_TYPE, 
          class FUNC_CLASS = defaultPileup<PILEUP_TYPE> >
class Pileup
{
public:
    Pileup(const FUNC_CLASS& fp = FUNC_CLASS());

    Pileup(int window, const FUNC_CLASS& fp = FUNC_CLASS());

    // Perform pileup with a reference.
    Pileup(const std::string& refSeqFileName,
           const FUNC_CLASS& fp = FUNC_CLASS());

    // Perform pileup with a reference.
    Pileup(int window, const std::string& refSeqFileName, 
           const FUNC_CLASS& fp = FUNC_CLASS());

    ~Pileup();

    // Performs a pileup on the specified file.
    // Returns 0 for success and non-zero for failure.
    // If excludeFlag is specified, any bit may be set for the
    // record to be dropped.
    // Defaulted to exclude:
    //    * unmapped, 
    //    * not primary alignment
    //    * failed platform/vendor quality check
    //    * PCR or optical duplicate
    // If includeFlag is specified, every bit must be set for the
    // record to be included - defaulted to 0 - no bits are required to be set.
    virtual int processFile(const std::string& fileName, 
                            uint16_t excludeFlag = 0x0704, 
                            uint16_t includeFlag = 0);

    // This should be overwritten to perform any necessary
    // filtering on the record.
    virtual void processAlignment(SamRecord& record);
   

    void flushPileup();

protected:
    FUNC_CLASS myAnalyzeFuncPtr;

    // Always need the reference position.
    void addAlignmentPosition(int refPosition, SamRecord& record);


    virtual void flushPileup(int refID, int refPosition);
    void flushPileup(int refPosition);
    
    // Get the position in the myElements container that is associated
    // with the specified position.  If the specified position cannot
    // fit within the myElements container, -1 is returned.
    int pileupPosition(int refPosition);

    virtual void resetElement(PILEUP_TYPE& element, int position);
    virtual void addElement(PILEUP_TYPE& element, SamRecord& record);
    virtual void analyzeElement(PILEUP_TYPE& element);
    virtual void analyzeHead();

    std::vector<PILEUP_TYPE> myElements;

    int    pileupStart;
    int    pileupHead;
    int    pileupTail;
    int    pileupWindow;

    int myCurrentRefID;

    GenomeSequence* myRefPtr;
};


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(1024),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    // Not using pointers since this is templated.
    myElements.resize(pileupWindow);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(int window, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(window),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    // Not using pointers since this is templated.
    myElements.resize(window);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(const std::string& refSeqFileName, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(1024),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    myRefPtr = new GenomeSequence(refSeqFileName.c_str());

    // Not using pointers since this is templated.
    myElements.resize(pileupWindow);

    PILEUP_TYPE::setReference(myRefPtr);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::Pileup(int window, const std::string& refSeqFileName, const FUNC_CLASS& fp)
    : myAnalyzeFuncPtr(fp),
      myElements(),
      pileupStart(0),
      pileupHead(0),
      pileupTail(-1),
      pileupWindow(window),
      myCurrentRefID(-2),
      myRefPtr(NULL)
{
    myRefPtr = new GenomeSequence(refSeqFileName.c_str());

    // Not using pointers since this is templated.
    myElements.resize(window);

    PILEUP_TYPE::setReference(myRefPtr);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
Pileup<PILEUP_TYPE, FUNC_CLASS>::~Pileup()
{
    if(myRefPtr != NULL)
    {
        delete myRefPtr;
        myRefPtr = NULL;
    }
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int Pileup<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& fileName, 
                                                 uint16_t excludeFlag,
                                                 uint16_t includeFlag)
{
    SamFile samIn;
    SamFileHeader header;
    SamRecord record;
    
    if(myRefPtr != NULL)
    {
        samIn.SetReference(myRefPtr);
    }

    if(!samIn.OpenForRead(fileName.c_str()))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }
    
    if(!samIn.ReadHeader(header))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

    // The file needs to be sorted by coordinate.
    samIn.setSortedValidation(SamFile::COORDINATE);

    // Iterate over all records
    while (samIn.ReadRecord(header, record))
    {
        uint16_t flag = record.getFlag();
        if(flag & excludeFlag)
        {
            // This record has an excluded flag set, 
            // so continue to the next one.
            continue;
        }
        if((flag & includeFlag) != includeFlag)
        {
            // This record does not have all required flags set, 
            // so continue to the next one.
            continue;
        }
        processAlignment(record);
    }

    flushPileup();

    int returnValue = 0;
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnValue = samIn.GetStatus();
    }
    return(returnValue);  
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::processAlignment(SamRecord& record)
{
    int refPosition = record.get0BasedPosition();
    int refID = record.getReferenceID();

    // Flush any elements from the pileup that are prior to this record
    // since the file is sorted, we are done with those positions.
    flushPileup(refID, refPosition);
    
    // Loop through for each reference position covered by the record.
    // It is up to the PILEUP_TYPE to handle insertions/deletions, etc
    // that are related with the given reference position.
    for(; refPosition <= record.get0BasedAlignmentEnd(); ++refPosition)
    {
        addAlignmentPosition(refPosition, record);
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup()
{
    // while there are still entries between the head and tail, flush,
    // but no need to flush if pileupTail == -1 because in that case 
    // no entries have been added
    while ((pileupHead <= pileupTail) && (pileupTail != -1))
    {
        flushPileup(pileupHead+1);
    }
    pileupStart = pileupHead = 0;
    pileupTail = -1;
}


// Always need the reference position.
template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::addAlignmentPosition(int refPosition,
                                                           SamRecord& record)
{
    int offset = pileupPosition(refPosition);
    
    if((offset < 0) || (offset >= pileupWindow))
    {
        std::cerr << "Pileup Buffer Overflow: position = " << refPosition
                  << "; refID = " << record.getReferenceID() 
                  << "; recStartPos = " << record.get1BasedPosition()
                  << "; pileupStart = " << pileupStart
                  << "; pileupHead = " << pileupHead
                  << "; pileupTail = " << pileupTail;
    }

    addElement(myElements[offset], record);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup(int refID, int position)
{
    // if new chromosome, flush the entire pileup.
    if(refID != myCurrentRefID)
    {
        // New chromosome, flush everything.
        flushPileup();
        myCurrentRefID = refID;
        // Set the starts to this position - first entry
        // in this chromosome.
        pileupStart = pileupHead = position;
    }
    else
    {
        // on the same chromosome, so flush just up to this new position.
        flushPileup(position);
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup(int position)
{
    // Flush up to this new position, but no reason to flush if
    // pileupHead has not been set.
    while((pileupHead < position) && (pileupHead <= pileupTail))
    {
        analyzeHead();

        pileupHead++;
        
        if(pileupHead - pileupStart >= pileupWindow)
            pileupStart += pileupWindow;
    }

    // If pileupHead != position, then we stopped because we had already
    // flushed all the data we had, just set the head and start to the position.
    if(pileupHead != position)
    {
        pileupHead = pileupStart = position;
    }
}


// Get the position in the myElements container that is associated
// with the specified position.  If the specified position cannot
// fit within the myElements container, -1 is returned.
template <class PILEUP_TYPE, class FUNC_CLASS>
int Pileup<PILEUP_TYPE, FUNC_CLASS>::pileupPosition(int position)
{
    if((position < pileupHead) || (position > (pileupHead + pileupWindow)))
    {
        return -1;
        error("Pileup Buffer Overflow");
    }    

    int offset = position - pileupStart;
    
    if(offset >= pileupWindow)
    {
        offset -= pileupWindow;
    }

    if(position > pileupTail)
    {
        // This is the first time this position is being used.
        // reset the element for this position.
        resetElement(myElements[offset], position);
        pileupTail = position;
    }

    return(offset);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::resetElement(PILEUP_TYPE& element,
                                                   int position)
{
    element.reset(position);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::addElement(PILEUP_TYPE& element,
                                                 SamRecord& record)
{
    element.addEntry(record);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::analyzeElement(PILEUP_TYPE& element)
{
    myAnalyzeFuncPtr(element);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void Pileup<PILEUP_TYPE, FUNC_CLASS>::analyzeHead()
{
    myAnalyzeFuncPtr(myElements[pileupHead - pileupStart]);
}


#endif
