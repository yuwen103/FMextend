///-----------------------------------------------
///-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ErrorCorrectProcess - Wrapper to perform error correction
// for a sequence work item
//
#include "ErrorCorrectProcess.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include "multiple_alignment.h"
#include "KmerOverlaps.h"
#include <iomanip>
#include "FMIndexWalkProcess.h"
#include "Extension.h"

//#define KMER_TESTING 1



//
//
//
ErrorCorrectProcess::ErrorCorrectProcess(const ErrorCorrectParameters params) : m_params(params)
{
	m_params.depthFilter = 10000;
}

//
ErrorCorrectProcess::~ErrorCorrectProcess()
{

}

ErrorCorrectResult ErrorCorrectProcess::process(const SequenceWorkItem& workItem)
{
        ErrorCorrectResult result = correct(workItem);
        if(!result.kmerQC && !result.overlapQC && m_params.printOverlaps)
        std::cout << workItem.read.id << " failed error correction QC\n";
        return result;
}

ErrorCorrectResult ErrorCorrectProcess::correct(const SequenceWorkItem& workItem)
{
	switch(m_params.algorithm)
	{
	case ECA_HYBRID:
		{
			ErrorCorrectResult result = kmerCorrection(workItem);
			if(!result.kmerQC)
			return overlapCorrectionNew(workItem);
			else
			return result;
			break;
		}
	case ECA_KMER:
		{
			return kmerCorrection(workItem);
			break;
		}
	case ECA_OVERLAP:
		{
			return overlapCorrectionNew(workItem);
			break;
		}
	case ECA_THREAD:
		{
			//return threadingCorrection(workItem);
			break;
		}
	default:
		{
			assert(false);
		}
	}
	ErrorCorrectResult result;
	return result;
}

//
ErrorCorrectResult ErrorCorrectProcess::overlapCorrectionNew(const SequenceWorkItem& workItem)
{
	assert(m_params.indices.pBWT != NULL);
	assert(m_params.indices.pSSA != NULL);

	ErrorCorrectResult result;
	std::string current_sequence = workItem.read.seq.toString();
	std::string consensus;
	
	
	
	//FE_tree_test
	std::string Query = current_sequence;
	
	
	int SolidKmer_size=m_params.kmerLength;
	size_t Seed_size2=m_params.check_kmerLength;
	//std::string read_id = workItem.read.id;
	//printf("The_ID=\t%s\n",read_id.c_str());
	//Extension::printfKFQ(Query,SolidKmer_size,m_params.indices.pBWT);
	//consensus=Extension::TrimReads(Query,SolidKmer_size,m_params.indices.pBWT);
	
	
	int k_diff=SolidKmer_size-(int)Seed_size2;
	Solid_error solid_info;
	std::vector<Ksub_vct> correct_ksub;
	std::vector<OutInfo> Out_Info;
	//printf("correct_ksub size=%d\n",(int)correct_ksub.size());
	correct_ksub.clear();
	//Extension::addQueryinKsub(Query,Seed_size2,correct_ksub);
	//Extension::addQueryinOutInfo(Query,Out_Info);
	/*
	std::string read_id = workItem.read.id;
	printf("The_ID=\t%s\n",read_id.c_str());
	Extension::printfKFQ(Query,SolidKmer_size,m_params.indices.pBWT);
	*/
	//printf("Query:\n%s\n",Query.c_str());
	
	solid_info=Extension::getSolidRegion(Query,SolidKmer_size,m_params.indices.pBWT);
	//printf("Final Start=\t%d\tEnd=\t%d\n",solid_info.solid_left_idx,solid_info.solid_right_idx);
	
	//printf("%d\t%d\t%s\t%s\n",solid_info.solid_left_idx,solid_info.solid_right_idx,solid_info.Left_error_1st_bp.c_str(),solid_info.Right_error_1st_bp.c_str());
	//consensus=Extension::getSolidRegion_v2(Query,SolidKmer_size,m_params.indices.pBWT);
	std::vector<BWTInterval> L_TerminatedIntervals,R_TerminatedIntervals; 
	L_TerminatedIntervals.clear();
	R_TerminatedIntervals.clear();
	std::string right_correct="";
	std::string left_correct="";
	std::string solid_Region="";
	std::string correct_str="";
	//printf("1 check\n");
	std::vector<Out_test> out_vct;
	out_vct.clear();
	
	if(solid_info.solid_left_idx>=0)
	{
		//std::string L_last_kmer=Query.substr(solid_info.solid_left_idx,(int)Seed_size2-1);
		//std::string R_last_kmer=Query.substr(solid_info.solid_right_idx+k_diff+1,(int)Seed_size2-1);
		
		//printf("2 check\n");
		Extension::getLRKmerInterval(Query,Seed_size2,m_params.indices.pBWT,m_params.indices.pRBWT,L_TerminatedIntervals,R_TerminatedIntervals);
		//printf("3 check\n");
		//printf("Left_size=%d\tRight_size=%d\n",(int)L_TerminatedIntervals.size(),(int)R_TerminatedIntervals.size());
		Extension::ExtensionRead(Query,SolidKmer_size,Seed_size2,solid_info.solid_right_idx,m_params.indices.pBWT,m_params.indices.pRBWT,Out_Info,L_TerminatedIntervals,R_TerminatedIntervals);
		//printf("Out_info_reads_count=\t%d\n",(int)Out_Info.size());
		//Extension::addStrInKsub(Out_Info,Seed_size2,correct_ksub);
		//Extension::addStrInKsub2(Out_Info,Seed_size2,out_vct,Query);
		
		/*
		for(int i=0;i<(int)Out_Info.size();i++)
		{
			std::string gap_temp="";
			for(int idx=0;idx<Out_Info[i].start_point;idx++)
			{
				gap_temp.append("=");
			}
			std::string out_str=gap_temp+Out_Info[i].outStr;
			printf(">%d\n%s\n",i,out_str.c_str());
			//printf("Str=\t%s\nSeed_rate=\t%.2f\tSize=\t%d\n",out_str.c_str(),Out_Info[i].Seed_rate,Out_Info[i].interval_size);
		}
		*/
		
		//int temp_kmer_size=15;
		
		
		//Seed_size2=1;
		//k_diff=SolidKmer_size-(int)Seed_size2;
		
		std::string L_last_kmer=Query.substr(solid_info.solid_left_idx,Seed_size2-1);
		std::string R_last_kmer=Query.substr(solid_info.solid_right_idx+(SolidKmer_size-Seed_size2)+1,Seed_size2-1);
		std::string L_First_kmer=Query.substr(solid_info.solid_left_idx,SolidKmer_size-1);
		std::string R_First_kmer=Query.substr(solid_info.solid_right_idx+1,SolidKmer_size-1);
		//test k-diff
		//k_diff=(SolidKmer_size-temp_kmer_size);
		
		//Extension::addStrInKsub3(Out_Info,Seed_size2,correct_ksub,Query);
		Extension::addStrInKsub3(Out_Info,Seed_size2,correct_ksub,Query);
		//printf("4 check\n");
		
		/*
		if(solid_info.solid_left_idx>0)
			left_correct=Extension::correct_left(solid_info,correct_ksub,L_last_kmer);
		right_correct=Extension::correct_right(solid_info,correct_ksub,k_diff,R_last_kmer);
		*/
		if(solid_info.solid_left_idx>0)
			left_correct=Extension::Newcorrect_left(Query,solid_info,correct_ksub,L_last_kmer,L_First_kmer,SolidKmer_size,m_params.indices.pBWT);
		right_correct=Extension::Newcorrect_right(Query,solid_info,correct_ksub,k_diff,R_last_kmer,R_First_kmer,SolidKmer_size,m_params.indices.pBWT);
		
		//printf("5 check\n");
		solid_Region=Query.substr(solid_info.solid_left_idx,(solid_info.solid_right_idx-solid_info.solid_left_idx+SolidKmer_size));
		correct_str=left_correct+solid_Region+right_correct;
		//printf("L=\t%s\nS=\t%s\nR=\t%s\n",left_correct.c_str(),solid_Region.c_str(),right_correct.c_str());
		//printf("6 check\n");
		
		//consensus=Extension::TrimReads(correct_str,31,m_params.indices.pBWT);
		
		//consensus=Extension::NewTrimReads(correct_str,31,m_params.indices.pBWT);
		
		//printf("Before:\t%s\n",correct_str.c_str());
		//printf("After:\t%s\n",consensus.c_str());
		
		//Extension::printfKFQ(correct_str,SolidKmer_size,m_params.indices.pBWT);
		//printf("Correct=\t%s\n",correct_str.c_str());
		//printf("Trim_st=\t%s\n",consensus.c_str());
		consensus=correct_str;
		
	}
	else if(solid_info.solid_left_idx==-3)
	{
		//printf("Didn't have any FQ>1 31-kmer.\n");
		consensus=current_sequence;
	}
	else
	{
		//printf("First round don't find Seed.\n");
		//int temp_solid_size=SolidKmer_size-10;
		int temp_solid_size=15;
		if(temp_solid_size>12)
		{
			Seed_size2=7;
			k_diff=temp_solid_size-(int)Seed_size2;
			
			//Extension::printfKFQ(Query,temp_solid_size,m_params.indices.pBWT);
			//solid_info=Extension::highError_getSolidRegion(Query,temp_solid_size,m_params.indices.pBWT);
			//Extension::printfKFQ(Query,temp_solid_size,m_params.indices.pBWT);
			solid_info=Extension::getSolidRegion(Query,temp_solid_size,m_params.indices.pBWT);
			
			if(solid_info.solid_left_idx>=0)
			{
				Extension::getLRKmerInterval(Query,Seed_size2,m_params.indices.pBWT,m_params.indices.pRBWT,L_TerminatedIntervals,R_TerminatedIntervals);
				Extension::ExtensionRead(Query,temp_solid_size,Seed_size2,solid_info.solid_right_idx,m_params.indices.pBWT,m_params.indices.pRBWT,Out_Info,L_TerminatedIntervals,R_TerminatedIntervals);
				//Seed_size2=1;
				//k_diff=SolidKmer_size-(int)Seed_size2;
				std::string L_last_kmer=Query.substr(solid_info.solid_left_idx,Seed_size2-1);
				std::string R_last_kmer=Query.substr(solid_info.solid_right_idx+(temp_solid_size-Seed_size2)+1,Seed_size2-1);
				std::string L_First_kmer=Query.substr(solid_info.solid_left_idx,temp_solid_size-1);
				std::string R_First_kmer=Query.substr(solid_info.solid_right_idx+1,temp_solid_size-1);
				
				//printf("The num of find reads = %d\n",(int)Out_Info.size());
				
				Extension::addStrInKsub3(Out_Info,Seed_size2,correct_ksub,Query);
				if(solid_info.solid_left_idx>0)
					left_correct=Extension::Newcorrect_left(Query,solid_info,correct_ksub,L_last_kmer,L_First_kmer,temp_solid_size,m_params.indices.pBWT);
				right_correct=Extension::Newcorrect_right(Query,solid_info,correct_ksub,k_diff,R_last_kmer,R_First_kmer,temp_solid_size,m_params.indices.pBWT);
				solid_Region=Query.substr(solid_info.solid_left_idx,(solid_info.solid_right_idx-solid_info.solid_left_idx+temp_solid_size));
				correct_str=left_correct+solid_Region+right_correct;
				//printf("Correct=\t%s\n",correct_str.c_str());
				//printf("Trim_st=\t%s\n",consensus.c_str());
				consensus=correct_str;
			}
		}
		else
		{
			//printf("No_correct\n");
		}
	}
	
	
	for(int test=0;test<(int)correct_ksub.size();test++)
	{
		printf("At position:%d Have Kmer:\n",test);
		for(int t=0;t<(int)correct_ksub[test].size();t++)
		{
			printf("%dth\t%s\tcount=\t%d\n",t,correct_ksub[test][t].kmer.c_str(),correct_ksub[test][t].countOfkmer);
		}
	}
	
	
	//printf("Left=\t%s\tSolid_region=\t%s\tRight=\t%s\n",left_correct.c_str(),
	//Query.substr(solid_info.solid_left_idx,(solid_info.solid_right_idx-solid_info.solid_left_idx+SolidKmer_size)).c_str(),right_correct.c_str());
	
	correct_ksub.clear();
	L_TerminatedIntervals.clear();
	R_TerminatedIntervals.clear();
	std::vector<Ksub_vct>().swap(correct_ksub);
	
	if(!consensus.empty())
	{
		result.correctSequence = consensus;
		result.overlapQC = true;
	}
	
	else
	{
		// Return the unmodified query sequence
		result.correctSequence = current_sequence;
		result.overlapQC = true;
	}
	
	return result;
}


// Correct a read with a k-mer based corrector
ErrorCorrectResult ErrorCorrectProcess::kmerCorrection(const SequenceWorkItem& workItem)
{
	assert(m_params.indices.pBWT != NULL);
	assert(m_params.indices.pCache != NULL);

	ErrorCorrectResult result;

	typedef std::map<std::string, int> KmerCountMap;
	KmerCountMap kmerCache;

	SeqRecord currRead = workItem.read;
	std::string readSequence = workItem.read.seq.toString();

#ifdef KMER_TESTING
	std::cout << "Kmer correcting read " << workItem.read.id << "\n";
#endif

	if((int)readSequence.size() < m_params.kmerLength)
	{
		// The read is shorter than the kmer length, nothing can be done
		result.correctSequence = readSequence;
		result.kmerQC = false;
		return result;
	}

	int n = readSequence.size();
	int nk = n - m_params.kmerLength + 1;

	// Are all kmers in the read well-represented?
	bool allSolid = false;
	bool done = false;
	int rounds = 0;
	int maxAttempts = m_params.numKmerRounds;

	// For each kmer, calculate the minimum phred score seen in the bases
	// of the kmer
	std::vector<int> minPhredVector(nk, 0);
	for(int i = 0; i < nk; ++i)
	{
		int end = i + m_params.kmerLength - 1;
		int minPhred = std::numeric_limits<int>::max();
		for(int j = i; j <= end; ++j)
		{
			int ps = workItem.read.getPhredScore(j);
			if(ps < minPhred)
			minPhred = ps;
		}
		minPhredVector[i] = minPhred;
	}

	while(!done && nk > 0)
	{
		// Compute the kmer counts across the read
		// and determine the positions in the read that are not covered by any solid kmers
		// These are the candidate incorrect bases
		std::vector<int> countVector(nk, 0);
		std::vector<int> solidVector(n, 0);

		for(int i = 0; i < nk; ++i)
		{
			std::string kmer = readSequence.substr(i, m_params.kmerLength);

			// First check if this kmer is in the cache
			// If its not, find its count from the fm-index and cache it
			int count = 0;
			KmerCountMap::iterator iter = kmerCache.find(kmer);

			if(iter != kmerCache.end())
			{
				count = iter->second;
			}
			else
			{
				count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
				kmerCache.insert(std::make_pair(kmer, count));
			}

			// Get the phred score for the last base of the kmer
			int phred = minPhredVector[i];
			countVector[i] = count;
			//            std::cout << i << "\t" << phred << "\t" << count << "\n";

			// Determine whether the base is solid or not based on phred scores
			int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);
			if(count >= threshold)
			{
				for(int j = i; j < i + m_params.kmerLength; ++j)
				solidVector[j] = 1;
			}
		}

		allSolid = true;
		for(int i = 0; i < n; ++i)
		{
#ifdef KMER_TESTING
			std::cout << "Position[" << i << "] = " << solidVector[i] << "\n";
#endif
			if(solidVector[i] != 1)
			allSolid = false;
		}

#ifdef KMER_TESTING
		std::cout << "Read " << workItem.read.id << (allSolid ? " is solid\n" : " has potential errors\n");
#endif

		// Stop if all kmers are well represented or we have exceeded the number of correction rounds
		if(allSolid || rounds++ > maxAttempts)
		break;

		// Attempt to correct the leftmost potentially incorrect base
		bool corrected = false;
		for(int i = 0; i < n; ++i)
		{
			if(solidVector[i] != 1)
			{
				// Attempt to correct the base using the leftmost covering kmer
				int phred = workItem.read.getPhredScore(i);
				int threshold = CorrectionThresholds::Instance().getRequiredSupport(phred);

				int left_k_idx = (i + 1 >= m_params.kmerLength ? i + 1 - m_params.kmerLength : 0);
				corrected = attemptKmerCorrection(i, left_k_idx, std::max(countVector[left_k_idx], threshold), readSequence);
				if(corrected)
				break;

				// base was not corrected, try using the rightmost covering kmer
				size_t right_k_idx = std::min(i, n - m_params.kmerLength);
				corrected = attemptKmerCorrection(i, right_k_idx, std::max(countVector[right_k_idx], threshold), readSequence);
				if(corrected)
				break;
			}
		}

		// If no base in the read was corrected, stop the correction process
		if(!corrected)
		{
			assert(!allSolid);
			done = true;
		}
	}

	if(allSolid)
	{
		result.correctSequence = readSequence;
		result.kmerQC = true;
	}
	else
	{
		result.correctSequence = workItem.read.seq.toString();
		result.kmerQC = false;
	}
	return result;
}


// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
// And there are other alleles with kmer freq <= avgCount and >= minCount
bool ErrorCorrectProcess::attemptHeteroCorrection(size_t i, size_t k_idx, size_t minCount, size_t avgCount, std::string& readSequence)
{
	assert(i >= k_idx && i < k_idx + m_params.kmerLength);
	size_t base_idx = i - k_idx;
	char originalBase = readSequence[i];

	std::string kmer = readSequence.substr(k_idx, m_params.kmerLength);
	int bestCount = -1;
	char bestBase = '$';

	bool isAnotherAlleleExisted=false;
	
	for(int j = 0; j < DNA_ALPHABET::size; ++j)
	{
		char currBase = ALPHABET[j];
		kmer[base_idx] = currBase;
		size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		//size_t count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);

		//Another allele must have kmer freq < avgCount and > minCount
		// std::cout << currBase << ":" << count << ":" << avgCount << "\n";
		if(count <= avgCount && count >= minCount*2)
		{
			if( currBase != originalBase)
				isAnotherAlleleExisted=true;
				
			if( (int)count > bestCount)
			{
				bestCount = count;
				bestBase = currBase;
			}
		}
	}

	if(isAnotherAlleleExisted)
	{
		readSequence[i] = bestBase;
		return true;
	}
	return false;
}

// Attempt to correct the base at position idx in readSequence. Returns true if a correction was made
// The correction is made only if the count of the corrected kmer is at least minCount
bool ErrorCorrectProcess::attemptKmerCorrection(size_t i, size_t k_idx, size_t minCount, std::string& readSequence)
{
	size_t kmerLength=m_params.kmerLength;
	assert(i >= k_idx && i < k_idx + kmerLength);
	size_t base_idx = i - k_idx;
	char originalBase = readSequence[i];
	std::string kmer = readSequence.substr(k_idx, kmerLength);
	size_t bestCount = 0;
	char bestBase = '$';

#if KMER_TESTING
	std::cout << "i: " << i << " k-idx: " << k_idx << " " << kmer << " " << reverseComplement(kmer) << "\n";
#endif

	for(int j = 0; j < DNA_ALPHABET::size; ++j)
	{
		char currBase = ALPHABET[j];
		// if(currBase == originalBase)
			// continue;
		kmer[base_idx] = currBase;
		size_t count = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		//size_t count = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);

// #if KMER_TESTING
		// printf("%c %c %zu\n", originalBase, currBase, count);
// #endif

		if(count >= minCount*2)
		{
			// printf("%c %c %zu\n", originalBase, currBase, count);
			// Multiple corrections exist, double kmer size
			// if(bestBase != '$' && kmerLength+m_params.kmerLength < readSequence.length())
			// {
				// kmerLength=kmerLength+m_params.kmerLength;
				// k_idx = (i >= kmerLength/2)?i - kmerLength/2:0;
				// if(k_idx>(int)readSequence.length()-kmerLength) k_idx=readSequence.length()-kmerLength;
				// base_idx = i - k_idx;
				// kmer = readSequence.substr(k_idx, kmerLength);	
				// bestCount = 0;
				// bestBase = '$';
				// j=-1;
				// continue;
			// }
			
			bestCount = count;
			bestBase = currBase;
		}
	}

	if(bestCount >= minCount*2 && bestBase!=originalBase)
	{
		assert(bestBase != '$');
		readSequence[i] = bestBase;
		return true;
	}
	return false;
}


//
//
//
ErrorCorrectPostProcess::ErrorCorrectPostProcess(std::ostream* pCorrectedWriter,
std::ostream* pDiscardWriter,
bool bCollectMetrics) :
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_bCollectMetrics(bCollectMetrics),
m_totalBases(0), m_totalErrors(0),
m_readsKept(0), m_readsDiscarded(0),
m_kmerQCPassed(0), m_overlapQCPassed(0),
m_qcFail(0)
{

}

//
ErrorCorrectPostProcess::~ErrorCorrectPostProcess()
{
	std::cout << "Reads passed kmer QC check: " << m_kmerQCPassed << "\n";
	std::cout << "Reads passed overlap QC check: " << m_overlapQCPassed << "\n";
	std::cout << "Reads failed QC: " << m_qcFail << "\n";
}

//
void ErrorCorrectPostProcess::writeMetrics(std::ostream* pWriter)
{
	m_positionMetrics.write(pWriter, "Bases corrected by position\n", "pos");
	m_originalBaseMetrics.write(pWriter, "\nOriginal base that was corrected\n", "base");
	m_precedingSeqMetrics.write(pWriter, "\nkmer preceding the corrected base\n", "kmer");
	m_qualityMetrics.write(pWriter, "\nBases corrected by quality value\n\n", "quality");

	std::cout << "ErrorCorrect -- Corrected " << m_totalErrors << " out of " << m_totalBases <<
	" bases (" << (double)m_totalErrors / m_totalBases << ")\n";
	std::cout << "Kept " << m_readsKept << " reads. Discarded " << m_readsDiscarded <<
	" reads (" << (double)m_readsDiscarded / (m_readsKept + m_readsDiscarded)<< ")\n";
}


//
void ErrorCorrectPostProcess::process(const SequenceWorkItem& item, const ErrorCorrectResult& result)
{

	// Determine if the read should be discarded
	bool readQCPass = true;
	if(result.kmerQC)
	{
		m_kmerQCPassed += 1;
	}
	else if(result.overlapQC)
	{
		m_overlapQCPassed += 1;
	}
	else
	{
		readQCPass = false;
		m_qcFail += 1;
	}


	// Collect metrics for the reads that were actually corrected
	if(m_bCollectMetrics && readQCPass)
	{
		collectMetrics(item.read.seq.toString(),
		result.correctSequence.toString(),
		item.read.qual);
	}

	SeqRecord record = item.read;
	record.seq = result.correctSequence;


	if (result.correctSequence.empty()) ;

	else if  (readQCPass || m_pDiscardWriter == NULL)
	{
		record.write(*m_pCorrectedWriter);
		++m_readsKept;
	}
	else
	{
		record.write(*m_pDiscardWriter);
		++m_readsDiscarded;
	}

}


void ErrorCorrectPostProcess::collectMetrics(const std::string& originalSeq,
const std::string& correctedSeq,
const std::string& qualityStr)
{
	size_t precedingLen = 2;
	for(size_t i = 0; i < originalSeq.length(); ++i)
	{
		char qc = !qualityStr.empty() ? qualityStr[i] : '\0';
		char ob = originalSeq[i];

		++m_totalBases;

		m_positionMetrics.incrementSample(i);

		if(!qualityStr.empty())
		m_qualityMetrics.incrementSample(qc);

		m_originalBaseMetrics.incrementSample(ob);

		std::string precedingMer;
		if(i > precedingLen)
		{
			precedingMer = originalSeq.substr(i - precedingLen, precedingLen);
			m_precedingSeqMetrics.incrementSample(precedingMer);
		}

		if(originalSeq[i] != correctedSeq[i])
		{
			m_positionMetrics.incrementError(i);
			if(!qualityStr.empty())
			m_qualityMetrics.incrementError(qc);
			m_originalBaseMetrics.incrementError(ob);

			if(!precedingMer.empty())
			{
				m_precedingSeqMetrics.incrementError(precedingMer);
			}
			++m_totalErrors;
		}
	}
}
