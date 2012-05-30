/* gfPcrLib - Routines to help do in silico PCR.
 * Copyright 2004-2005 Jim Kent.  All rights reserved. */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dystring.h"
#include "fa.h"
#include "net.h"
#include "genoFind.h"
#include "sqlNum.h"
#include "gfInternal.h"
#include "gfPcrLib.h"

static char const rcsid[] = "$Id: gfPcrLib.c,v 1.8 2006/03/15 18:36:16 angie Exp $";

/**** Input and Output Handlers *****/

void gfPcrOutputFree(struct gfPcrOutput **pOut)
/* Free up a gfPcrOutput structure. */
{
    struct gfPcrOutput *out = *pOut;
    if (pOut != NULL)
    {
        freeMem(out->name);
        freeMem(out->fPrimer);
        freeMem(out->rPrimer);
        freeMem(out->seqName);
        freeMem(out->dna);
        freez(pOut);
    }
}

void gfPcrOutputFreeList(struct gfPcrOutput **pList)
/* Free up a list of gfPcrOutput structures. */
{
    struct gfPcrOutput *el, *next;

    for (el = *pList; el != NULL; el = next)
    {
        next = el->next;
        gfPcrOutputFree(&el);
    }
    *pList = NULL;
}

void gfPcrInputStaticLoad(char **row, struct gfPcrInput *ret)
/* Load a row from gfPcrInput table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{
    ret->name = row[0];
    ret->fPrimer = row[1];
    ret->rPrimer = row[2];
}

struct gfPcrInput *gfPcrInputLoad(char **row)
/* Load a gfPcrInput from row fetched with select * from gfPcrInput
 * from database.  Dispose of this with gfPcrInputFree(). */
{
    struct gfPcrInput *ret;

    AllocVar(ret);
    ret->name = cloneString(row[0]);
    ret->fPrimer = cloneString(row[1]);
    ret->rPrimer = cloneString(row[2]);
    return ret;
}

struct gfPcrInput *gfPcrInputLoadAll(char *fileName)
/* Load all gfPcrInput from a whitespace-separated file.
 * Dispose of this with gfPcrInputFreeList(). */
{
    struct gfPcrInput *list = NULL, *el;
    struct lineFile *lf = lineFileOpen(fileName, TRUE);
    char *row[3];

    while (lineFileRow(lf, row))
    {
        el = gfPcrInputLoad(row);
        slAddHead(&list, el);
    }
    lineFileClose(&lf);
    slReverse(&list);
    return list;
}


void gfPcrInputFree(struct gfPcrInput **pEl)
/* Free a single dynamically allocated gfPcrInput such as created
 * with gfPcrInputLoad(). */
{
    struct gfPcrInput *el;

    if ((el = *pEl) == NULL) return;
    freeMem(el->name);
    freeMem(el->fPrimer);
    freeMem(el->rPrimer);
    freez(pEl);
}

void gfPcrInputFreeList(struct gfPcrInput **pList)
/* Free a list of dynamically allocated gfPcrInput's */
{
    struct gfPcrInput *el, *next;

    for (el = *pList; el != NULL; el = next)
    {
        next = el->next;
        gfPcrInputFree(&el);
    }
    *pList = NULL;
}



/**** Handle Refinement (after Index has given us approximate position) *****/

static boolean goodMatch(char *a, char *b, int size)
/* Return TRUE if there are 2 matches per mismatch. */
{
    int score = 0, i;
    for (i=0; i<size; ++i)
    {
        if (a[i] == b[i])
            score += 1;
        else
            score -= 2;
    }
    return score >= 0;
}

static void upperMatch(char *dna, char *primer, int size)
/* Uppercase DNA where it matches primer. */
{
    int i;
    for (i=0; i<size; ++i)
    {
        if (dna[i] == primer[i])
            dna[i] = toupper(dna[i]);
    }
}

static void outputFa(struct gfPcrOutput *out, FILE *f, char *url)
/* Output match in fasta format. */
{
    int fPrimerSize = strlen(out->fPrimer);
    int rPrimerSize = strlen(out->rPrimer);
    int productSize = out->rPos - out->fPos;
    char *dna = cloneStringZ(out->dna, productSize);
    char *rrPrimer = cloneString(out->rPrimer);
    char *ffPrimer = cloneString(out->fPrimer);
    struct dyString *faLabel = newDyString(0);
    char *name = out->name;

    /* Create fasta header with position, possibly empty name, and upper cased primers with position optionally hyperlinked. */
    touppers(rrPrimer);
    touppers(ffPrimer);
    if (url != NULL)
    {
        dyStringAppend(faLabel, "<A HREF=\"");
        dyStringPrintf(faLabel, url, out->seqName, out->fPos+1, out->rPos);
        dyStringAppend(faLabel, "\">");
    }
    dyStringPrintf(faLabel, "%s:%d%c%d",
                   out->seqName, out->fPos+1, out->strand, out->rPos);
    if (url != NULL)
        dyStringAppend(faLabel, "</A>");
    if (name != NULL)
        dyStringPrintf(faLabel, " %s", name);
    dyStringPrintf(faLabel, " %dbp %s %s",
                   out->rPos - out->fPos, ffPrimer, rrPrimer);

    /* Flip reverse primer to be in same direction and case as sequence. */
    reverseComplement(rrPrimer, rPrimerSize);
    tolowers(rrPrimer);

    /* Capitalize where sequence and primer match, and write out sequence. */
    upperMatch(dna, out->fPrimer, fPrimerSize);
    upperMatch(dna + productSize - rPrimerSize, rrPrimer, rPrimerSize);
    faWriteNext(f, faLabel->string, dna, productSize);

    /* Clean up. */
    freez(&dna);
    freez(&rrPrimer);
    freez(&ffPrimer);
    dyStringFree(&faLabel)
}

static int countMatch(char *a, char *b, int size)
/* Count number of letters in a, b that match */
{
    int count = 0, i;
    for (i=0; i<size; ++i)
        if (a[i] == b[i])
            ++count;
    return count;
}

static void outputBed(struct gfPcrOutput *out, FILE *f, char *url)
/* Output match in BED format. */
{
    int match;
    int size = out->rPos - out->fPos;
    int fPrimerSize = strlen(out->fPrimer);
    int rPrimerSize = strlen(out->rPrimer);
    char *name = out->name;
    if (name == NULL) name = "n/a";
    fprintf(f, "%s\t%d\t%d\t", out->seqName, out->fPos, out->rPos);
    fprintf(f, "%s\t", name);
    match = countMatch(out->dna, out->fPrimer, fPrimerSize);
    assert(size > 0);
    reverseComplement(out->rPrimer, rPrimerSize);
    match += countMatch(out->dna + size - rPrimerSize, out->rPrimer, rPrimerSize);
    reverseComplement(out->rPrimer, rPrimerSize);
    fprintf(f, "%d\t", round(1000.0 * match / (double)(fPrimerSize + rPrimerSize) ));
    fprintf(f, "%c\n", out->strand);
}

static void outputPsl(struct gfPcrOutput *out, FILE *f, char *url)
/* Output match in PSL format. */
{
    int match;
    int size = out->rPos - out->fPos;
    int fPrimerSize = strlen(out->fPrimer);
    int rPrimerSize = strlen(out->rPrimer);
    int bothSize = fPrimerSize + rPrimerSize;
    int gapSize = size - bothSize;
    char *name = out->name;
    if (name == NULL) name = "n/a";
    match = countMatch(out->dna, out->fPrimer, fPrimerSize);
    reverseComplement(out->rPrimer, rPrimerSize);
    assert(size > 0);
    match += countMatch(out->dna + size - rPrimerSize, out->rPrimer, rPrimerSize);
    reverseComplement(out->rPrimer, rPrimerSize);

    fprintf(f, "%d\t", match);
    fprintf(f, "%d\t", bothSize - match);
    fprintf(f, "0\t0\t");	/* repMatch, nCount. */
    fprintf(f, "1\t%d\t", gapSize);   /* qNumInsert, qBaseInsert */
    fprintf(f, "1\t%d\t", gapSize);   /* tNumInsert, tBaseInsert */
    fprintf(f, "%c\t", out->strand);
    fprintf(f, "%s\t", name);
    fprintf(f, "%d\t", size);
    fprintf(f, "0\t%d\t", size);	/* qStart, qEnd */
    fprintf(f, "%s\t%d\t", out->seqName, out->seqSize);
    fprintf(f, "%d\t%d\t", out->fPos, out->rPos);
    fprintf(f, "2\t");
    fprintf(f, "%d,%d,\t", fPrimerSize, rPrimerSize);
    fprintf(f, "%d,%d,\t", 0,size - rPrimerSize);
    fprintf(f, "%d,%d,\n", out->fPos, out->rPos - rPrimerSize);
}


void gfPcrOutputWriteList(struct gfPcrOutput *outList, char *outType,
                          char *url, FILE *f)
/* Write list of outputs in specified format (either "fa" or "bed")
 * to file.  If url is non-null it should be a printf formatted
 * string that takes %s, %d, %d for chromosome, start, end. */
{
    struct gfPcrOutput *out;
    void (*output)(struct gfPcrOutput *out, FILE *f, char *url) = NULL;
    if (sameWord(outType, "fa"))
        output = outputFa;
    else if (sameWord(outType, "bed"))
        output = outputBed;
    else if (sameWord(outType, "psl"))
        output = outputPsl;
    else
        errAbort("Unrecognized pcr output type %s", outType);
    for (out = outList; out != NULL; out = out->next)
        output(out, f, url);
}

void gfPcrOutputWriteAll(struct gfPcrOutput *outList,
                         char *outType, char *url, char *fileName)
/* Create file of outputs in specified format (either "fa" or "bed")
 * to file.  If url is non-null it should be a printf formatted
 * string that takes %s, %d, %d for chromosome, start, end. */
{
    FILE *f = mustOpen(fileName, "w");
    gfPcrOutputWriteList(outList, outType, url, f);
    carefulClose(&f);
}


static void pcrLocalStrand(char *pcrName,
                           struct dnaSeq *seq,  int seqOffset, char *seqName, int seqSize,
                           int maxSize, char *fPrimer, int fPrimerSize, char *rPrimer, int rPrimerSize,
                           int minPerfect, int minGood,
                           char strand, struct gfPcrOutput **pOutList)
/* Do detailed PCR scan on one strand and report results. */
{
    char *fDna = seq->dna, *rDna;
    char *endDna = fDna + seq->size;
    char *fpPerfect = fPrimer + fPrimerSize - minPerfect;
    char *rpPerfect = rPrimer;
    char *fpMatch, *rpMatch;
    int sizeLeft;
    int goodSize = minGood - minPerfect;
    struct gfPcrOutput  *out;
    int matchSize;

    reverseComplement(rPrimer, rPrimerSize);
    for (;;)
    {
        fpMatch = memMatch(fpPerfect, minPerfect, fDna, endDna - fDna);
        if (fpMatch == NULL)
            break;
        sizeLeft = endDna - fpMatch;
        rDna = fpMatch;
        for (;;)
        {
            int fPos, rPos, fGoodPos, rGoodPos;
            rpMatch = memMatch(rpPerfect, minPerfect, rDna, endDna - rDna);
            if (rpMatch == NULL)
                break;
            fPos = fpMatch - (fPrimerSize - minPerfect) - seq->dna;
            rPos = rpMatch + rPrimerSize - seq->dna;
            fGoodPos = fpMatch - goodSize - seq->dna;
            rGoodPos = rpMatch + minPerfect - seq->dna;
            matchSize = rPos - fPos;
            if (rPos >= 0 && fPos >= 0 && fPos < seq->size && matchSize <= maxSize)
            {
                /* If matches well enough create output record. */
                if (goodMatch(seq->dna + fGoodPos, fpPerfect - goodSize, goodSize) &&
                        goodMatch(seq->dna + rGoodPos, rpPerfect + minPerfect, goodSize))
                {
                    AllocVar(out);
                    out->name  = cloneString(pcrName);
                    out->fPrimer = cloneString(fPrimer);
                    out->rPrimer = cloneString(rPrimer);
                    reverseComplement(out->rPrimer, rPrimerSize);
                    out->seqName = cloneString(seqName);
                    out->seqSize = seqSize;
                    out->fPos = fPos + seqOffset;
                    out->rPos = rPos + seqOffset;
                    out->strand = strand;
                    out->dna = cloneStringZ(seq->dna + fPos, matchSize);

                    /* Dealing with the strand of darkness....  Here we just have to swap
                     * forward and reverse primers to flip strands, and reverse complement
                     * the amplified area.. */
                    if (strand == '-')
                    {
                        char *temp = out->rPrimer;
                        out->rPrimer = out->fPrimer;
                        out->fPrimer = temp;
                        reverseComplement(out->dna, matchSize);
                    }
                    slAddHead(pOutList, out);
                }
            }
            rDna = rpMatch+1;
        }
        fDna = fpMatch + 1;
    }
    reverseComplement(rPrimer, rPrimerSize);
}

void gfPcrLocal(char *pcrName,
                struct dnaSeq *seq, int seqOffset, char *seqName, int seqSize,
                int maxSize, char *fPrimer, int fPrimerSize, char *rPrimer, int rPrimerSize,
                int minPerfect, int minGood, char strand, struct gfPcrOutput **pOutList)
/* Do detailed PCR scan on DNA already loaded into memory and put results
 * (in reverse order) on *pOutList. */
{
    /* For PCR primers reversing search strand just means switching
     * order of primers. */
    if (strand == '-')
        pcrLocalStrand(pcrName, seq, seqOffset, seqName, seqSize, maxSize,
                       rPrimer, rPrimerSize, fPrimer, fPrimerSize,
                       minPerfect, minGood, strand, pOutList);
    else
        pcrLocalStrand(pcrName, seq, seqOffset, seqName, seqSize, maxSize,
                       fPrimer, fPrimerSize, rPrimer, rPrimerSize,
                       minPerfect, minGood, strand, pOutList);
}

struct gfRange *gfPcrGetRanges(char *host, char *port, char *fPrimer, char *rPrimer,
                               int maxSize)
/* Query gfServer with primers and convert response to a list of gfRanges. */
{
    char buf[256];
    int conn = gfConnect(host, port);
    struct gfRange *rangeList = NULL, *range;

    /* Query server and put results into rangeList. */
    safef(buf, sizeof(buf), "%spcr %s %s %d", gfSignature(), fPrimer, rPrimer, maxSize);
    write(conn, buf, strlen(buf));
    for (;;)
    {
        if (netGetString(conn, buf) == NULL)
            break;
        if (sameString(buf, "end"))
            break;
        else if (startsWith("Error:", buf))
            errAbort(buf);
        else
        {
            char *s = buf;
            char *name, *start, *end, *strand;
            name = nextWord(&s);
            start = nextWord(&s);
            end = nextWord(&s);
            strand = nextWord(&s);
            if (strand == NULL)
                errAbort("Truncated gfServer response");
            AllocVar(range);
            range->tName = cloneString(name);
            range->tStart = atoi(start);
            range->tEnd = atoi(end);
            range->tStrand = strand[0];
            slAddHead(&rangeList, range);
        }
    }
    close(conn);
    slReverse(&rangeList);
    return rangeList;
}

static void gfPcrOneViaNet(
    char *host, char *port, char *seqDir,
    char *pcrName, char *fPrimer, char *rPrimer, int maxSize,
    int minPerfect, int minGood,
    struct hash *tFileCache, struct gfPcrOutput **pOutList)
/* Query gfServer for likely primer hits, load DNA to do detailed
 * examination, and output hits to head of *pOutList.. */
{
    struct gfRange *rangeList = NULL, *range;
    int fPrimerSize = strlen(fPrimer);
    int rPrimerSize = strlen(rPrimer);
    int maxPrimerSize = max(fPrimerSize, rPrimerSize);
    int minPrimerSize = min(fPrimerSize, rPrimerSize);

    tolowers(fPrimer);
    tolowers(rPrimer);

    /* Check for obvious user mistake. */
    if (minPrimerSize < minGood || minPrimerSize < minPerfect)
        errAbort("Need longer primer in pair %s (%dbp) %s (%dbp).  Minimum size is %d\n",
                 fPrimer, fPrimerSize, rPrimer, rPrimerSize, minGood);

    /* Load ranges and do more detailed snooping. */
    rangeList = gfPcrGetRanges(host, port, fPrimer, rPrimer, maxSize);
    for (range = rangeList; range != NULL; range = range->next)
    {
        int tSeqSize;
        char seqName[PATH_LEN];
        struct dnaSeq *seq = gfiExpandAndLoadCached(range,
                             tFileCache, seqDir,  0, &tSeqSize, FALSE, FALSE, maxPrimerSize);
        gfiGetSeqName(range->tName, seqName, NULL);
        gfPcrLocal(pcrName, seq, range->tStart, seqName, tSeqSize, maxSize,
                   fPrimer, fPrimerSize, rPrimer, rPrimerSize,
                   minPerfect, minGood, range->tStrand, pOutList);
        dnaSeqFree(&seq);
    }
    gfRangeFreeList(&rangeList);
}


struct gfPcrOutput *gfPcrViaNet(char *host, char *port, char *seqDir,
                                struct gfPcrInput *inList, int maxSize, int minPerfect, int minGood)
/* Do PCRs using gfServer index, returning list of results. */
{
    struct hash *tFileCache = gfFileCacheNew();
    struct gfPcrInput *in;
    struct gfPcrOutput *outList = NULL;

    for (in = inList; in != NULL; in = in->next)
    {
        gfPcrOneViaNet(host, port, seqDir,
                       in->name, in->fPrimer, in->rPrimer, maxSize,
                       minPerfect, minGood,
                       tFileCache, &outList);
    }
    gfFileCacheFree(&tFileCache);
    slReverse(&outList);
    return outList;
}

char *gfPcrMakePrimer(char *s)
/* Make primer (lowercased DNA) out of text.  Complain if
 * it is too short or too long. */
{
    int size = dnaFilteredSize(s);
    int realSize;
    char *primer = needMem(size+1);
    dnaFilter(s, primer);
    realSize = size - countChars(primer, 'n');
    if (realSize < 10 || realSize < size/2)
        errAbort("%s does not seem to be a good primer", s);
    return primer;
}

