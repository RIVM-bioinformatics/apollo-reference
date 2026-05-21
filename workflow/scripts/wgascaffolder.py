from __future__ import annotations

# Python Imports
import sys
import os
import pysam
from dataclasses import dataclass, field
from typing import Optional, Any, Union, List
from types import SimpleNamespace
from pathlib import Path

# Python Imports; needed for reporting
import portion as P
from tabulate import tabulate
import pandas as pd

# --- CORE DATA STRUCTURES ---

class SegmentArithmeticMixin:
    """ """
    def _get_overlap(self, s2:AnySegmentType) -> int:
        """ calculate overlap between two GenericSegments """
        s1 = self
        if not s1 or not s2: return 0
        if s1.target != s2.target: return 0
        #return max(0, min(s1[1], s2[1]) - max(s1[0], s2[0]))
        return max(0, min(s1.end, s2.end) - max(s1.start, s2.start))

    def _get_shortest_target_length(self) -> int:
        """ get shortest involved segment length for arithmetic on relative length(s) """
        if hasattr(self,'target_len'):
            return self.target_len
        else:
            lengths = list(sorted([self.ref.target_len,self.query.target_len]))
            if lengths[0] in (0,None):
                return lengths[1]
            else:
                return lengths[0]

    def _get_absolute_overlap(self,s2:GenericSegment) -> int:
        return self._get_overlap(s2)

    def _get_relative_overlap(self, s2: GenericSegment) -> float:
        ovl = self._get_overlap(s2)
        _length = self._get_shortest_target_length()
        #print(ovl,self.length,s2.length,_length)
        return float(ovl)/_length

@dataclass
class GenericSegment(SegmentArithmeticMixin):
    target: str
    start: int
    end: int
    name: str
    strand: str
    target_len: int
    length: int = field(init=False)

    def __post_init__(self):
        self.length = self.end - self.start

@dataclass
class AlignedSegment(SegmentArithmeticMixin):
    ref: GenericSegment
    query: GenericSegment
    score: float
    idx: int  # 1-based index
    raw: Any
    #phase: int = -1  # 0...N: Main, -2: Nested, -1: Bumped
    length: int = field(init=False)
    target: str = field(init=False)
    start: int = field(init=False)
    end: int = field(init=False)
    name: str = field(init=False)

    def __post_init__(self):
        # self.length != alignment_length.
        # but aligned fraction of reference and query can be distinct in length
        # here it's best to take the longest length, since it will be closer to the alignment_length
        self.length = max([ self.ref.end - self.ref.start, self.query.end - self.query.start ])
        self.target = self.ref.target
        self.start = self.ref.start
        self.end = self.ref.end
        self.name = self.ref.name


# typehint
AnySegmentType = Optional[Union[GenericSegment, AlignedSegment]]

class InputDataProcessor:
    def _expose_inputdata(self, file_or_filepointer:Any) -> Any:
        """ """
        msg = "please define in YourMatrixScaffolder(MatrixScaffolder) an implemented _expose_inputdata() class method"
        raise NotImplementedError(msg)

    def _loop_over_inputdata(self) -> None:
        """ """
        msg = "please define in YourMatrixScaffolder(MatrixScaffolder) an implemented _loop_over_inputdata() class method"
        raise NotImplementedError(msg)

    def _wrap(self, item: Any, idx: int) -> AnySegmentType:
        """ """
        msg = "please define in YourMatrixScaffolder(MatrixScaffolder) an implemented _wrap() class method"
        raise NotImplementedError(msg)

class OutputDataProcessor:
    def _prepare_output_files(self):
        """ hook to (e.g.) open various/complex output files """
        # optionally override. Default behaviour is doing nothing
        pass

    def _flush_scaffolded_input_data(self, out=Any):
        """ flush the (raw, putative memory-heavy) input data into your custom desired output """
        # optionally override. Default behaviour is doing nothing
        pass

    def _finalize_output_files(self):
        """ hook to (e.g.) close or stdout.write various/complex output files """
        # optionally override. Default behaviour is doing nothing
        pass

    def configure_output(self, **kwargs):
        """ provide configuration into the instance which the other OutputDataProcessor class methods might require """
        self.output_config.update(kwargs)
        return self

# MatrixScaffolder defaults (needed easier for editing configuration)
defaults = dict(
    ploidy=1,
    min_nt_length=None,
    min_unique_contribution=10000,
    max_overlap_abs=2500,
    max_overlap_rel=0.1
    )

MatrixScaffolderDefaults = SimpleNamespace(**defaults)
MSD = MatrixScaffolderDefaults

class MatrixScaffolder(InputDataProcessor,OutputDataProcessor):
    def __init__(self,
        ploidy:int=MSD.ploidy,
        min_nt_length:int=MSD.min_nt_length,
        min_unique_contribution:int=MSD.min_unique_contribution,
        max_overlap_abs:int=MSD.max_overlap_abs,
        max_overlap_rel:float=MSD.max_overlap_rel
    ):
        self.ploidy = ploidy
        self.max_overlap_abs = max_overlap_abs
        self.max_overlap_rel = max_overlap_rel
        self.min_unique_contribution = min_unique_contribution
        self.min_nt_length = min_nt_length

        # The MatrixScaffold attribute itself
        self.active_phases:tuple = tuple()
        self._reset_empty_active_phases()

        # stored information when flushing active_phases upon new contig
        self.data_summary = []
        # placeholder for complex output generation. Typically not used / overkill
        self.output_config = {}

        # might/will be needed for future cases
        #self.nested_pool: List[AlignedSegment] = []
        #self.skipped_pool: List[AlignedSegment] = []

        # current segment pointers
        self.current_chrom:str = None
        self.last_pos:int = 0

    @property
    def current_segment(self) -> AnySegmentType:
        try:
            return self.active_phases[0][-1]
        except IndexError:
            return None

    def get_segments_by_rowId(self,rowId:int) -> List[AnySegmentType]:
        """ """
        segments = []
        for phase in range(self.ploidy):
            segm = self.active_phases[phase][rowId]
            if not segm:
                continue
            segments.append(segm)
        return segments

    def run(self, inputdata:Any):
        """ public API to process inputdata into outputdata with all the IO hooks """
        self.inputdata = self._expose_inputdata(inputdata)
        self._prepare_output_files()
        try:
            self._loop_over_inputdata()
        finally:
            # !important! process final contig
            self._flush_current_chrom_phase_table()
        self._finalize_output_files()

    def _reset_empty_active_phases(self) -> None:
        """ """
        self.active_phases = tuple([list() for _ in range(self.ploidy + 1)])

    def has_vacant_phase(self) -> bool:
        """ """
        if self.current_chrom == None:
            return True
        elif self.ploidy == 1:
            return False
        elif True in [ self.active_phases[phase][-1] == None for phase in range(self.ploidy) ]:
            return True
        else:
            return False

    def _post_process_segment(self,segm:AnySegmentType) -> None:
        """ """
        self.current_chrom = segm.ref.target
        self.last_pos = max([self.last_pos,segm.ref.end])

    def _place_segment_on_new_row_hook(self,segm:AnySegmentType) -> None:
        """ """
        for i, ph in enumerate(self.active_phases):
            self.active_phases[i].append(None)
        # final "phase" allows multiple segments!
        self.active_phases[-1][-1] = []
        self.active_phases[0][-1] = segm

    def _place_segment_on_existing_row_hook(self,segm:AnySegmentType) -> None:
        """ """
        ## get current segments in one row
        #current_segments = [ self.active_phases[idx] for idx in range(len(self.active_phases)) ]
        #current_segments = current_segments[0:-1] + current_segments[-1]
        #while None in current_segments:
        #    current_segments.pop(current_segments.index(None))

        exceeds_length = None
        for i, phase in enumerate(self.active_phases):
            if phase[-1] == None:
                # vacant position!
                phase[-1] = segm
                if self.ploidy >= 2 and i >= 1:
                    # rule of thumb: try to populate the lowest phase with the largest segment
                    segments = [ self.active_phases[j][-1] for j in range(0,i+1) ]
                    segments.sort(key=lambda x: x.length, reverse=True)
                    for j,seg in enumerate(segments):
                        self.active_phases[j][-1] = seg
                # segment is placed. for .. else statement will not be activated!
                break
            elif i == len(self.active_phases)-1:
                # phase "bumped" is a list (of segments), not a segment
                pass
            else:
                exceeds_length = segm.length > phase[-1].length
        else:
            # no phase vacant? check if segment length exceeds earlier placed segments
            if exceeds_length:
                for i, phase in enumerate(self.active_phases[0:-1]):
                    if segm.length > phase[-1].length:
                        # Current segment is preferred over an earlier assigned segment.
                        # Bump the earlier one from its chair and take its place
                        # Rule of thumb: try to populate the lowest phase with the largest segment
                        # Since all chairs are taken, after segment sorting,
                        # the final segment (shortest/lowestpriority) will get bumped out of all chairs
                        segments = [self.active_phases[j][-1] for j in range(0, self.ploidy)]
                        segments.append(segm)
                        segments.sort(key=lambda x: x.length, reverse=True)
                        bumped = segments.pop()
                        for j, seg in enumerate(segments):
                            self.active_phases[j][-1] = seg
                        self.active_phases[-1][-1].append(bumped)
                        break
                else:
                    # append to the "bumped" phase
                    self.active_phases[-1][-1].append(segm)
            else:
                self.active_phases[-1][-1].append(segm)

            # Inspect self.active_phases[-1][-1], being the last bumped segment
            # Evaluate if it should be rescued because of min_unique_contribution
            # Realize, based on above logics, the "bumped" segment can be:
            # - current segment (the provided variable)
            # - or the one which was freshly bumped by current segment!
            bumped = self.active_phases[-1][-1][-1]

            # Although overlapping (otherwise we'd not have ended up here),
            # check if segment spans a large enough unique area
            # In case this is true, we have to overrule _place_segment_on_existing_row_hook
            # Fallback to _place_segment_on_new_row_hook
            if not self.min_unique_contribution:
                pass
            elif bumped.end - self.last_pos > self.min_unique_contribution:
                self.active_phases[-1][-1].pop()
                self._place_segment_on_new_row_hook(bumped)
            elif (self.current_segment.start - bumped.start) >= self.min_unique_contribution:
                ####print("##",segm.idx,bumped.idx,[ s.idx for s in self.active_phases[-1][-1]])
                self.active_phases[-1][-1].pop()
                self._place_segment_on_new_row_hook(bumped)
                # move upwards one element in the queue (bumped is on front of self.current_segment)
                for phase in self.active_phases:
                    _segm = phase.pop()
                    phase.insert(-1,_segm)
                # don't forget to move the "bumped" ones of this element along one position upwards
                self.active_phases[-1][-2] = self.active_phases[-1][-1]
                self.active_phases[-1][-1] = []
            else:
                # failing min_unique_contribution threshold
                pass

    def _flush_current_chrom_phase_table(self):
        """ Flush current chromosome: augment summary to self.data_summary and _flush_input_data """

        data = []
        for rowId in range(len(self.active_phases[0])):
            if self.active_phases[0][rowId].target != self.current_chrom:
                continue
            assigned_area = P.empty()
            bumped_area = P.empty()
            bumped_len = 0
            any_bumped_segment = False
            segm = self.active_phases[0][rowId]
            row = [segm.target,segm.start,segm.end,segm.name,segm.length,"+",'.']
            if len(data) > 0:
                segm = self.active_phases[0][rowId]
                dist = {True:"."}.get(segm.target != data[-1][0],segm.start - data[-1][2])
                row[-1] = dist

            for phase in range(self.ploidy+1):
                segm = self.active_phases[phase][rowId]
                if type(segm) == list and len(segm) == 0:
                    row.append(".")
                elif type(segm) == list and len(segm) >= 1:
                    any_bumped_segment = True
                    row.append(",".join([ str(s.idx) for s in segm ]) or ".")
                    for s in segm:
                        bumped_area |= P.closed(s.start, s.end)
                        bumped_len += ( s.end - s.start )
                elif segm == None:
                    row.append(".")
                else:
                    # one of the main phases
                    row.append(segm.idx)
                    assigned_area |= P.closed(segm.start, segm.end)

            if self.ploidy >= 2:
                row.insert(-1,".")
                if self.active_phases[1][rowId]:
                    row[-2] = self.active_phases[1][rowId].name
                segments = self.get_segments_by_rowId(rowId)
                _min = min([ segm.start for segm in segments ])
                _max = max([ segm.end for segm in segments ])
                row = [ segments[0].target, _min, _max ] + row


            if any_bumped_segment:
                # check size of the bumped area
                dropped_area = bumped_area - assigned_area
                subranges = [(interval.lower, interval.upper) for interval in dropped_area]
                length_uniq = sum([u-l for (l,u) in subranges])
                #print("XX",subranges,length,len(dropped_area))
                row.insert(-1,length_uniq)
                row.insert(-1,bumped_len)
                if len(row[-1]) >= 10:
                    row[-1] = "n=%s" % len(row[-1])
            else:
                row.insert(-1,0)
            data.append(row)

        # Reported. Now active_phases "queue" can be cleaned/reset
        self._flush_scaffolded_input_data()
        self.data_summary.extend(data)
        self._reset_empty_active_phases()

    def report(self,format:str='fancy') -> None:    # -> pd.DataFrame:
        if format == 'tsv':
            for row in self.data_summary:
                print("\t".join(map(str,row)))
        elif format == 'fancy':
            header = "chr start end name length strand dist".split()
            # TODO: specify oxaploid is the highest supported ploidy
            for phase in ("I","II","III","IV","V","VI","VII","VIII")[0:self.ploidy]:
                header.append("phase_"+phase)
            if self.ploidy >= 2:
                header.append("other")
                header = ["CHR","STA","END"] + header
            header.extend(['dropped_nt','dropped_len','dropped_idx'])
            df = pd.DataFrame.from_records(self.data_summary, columns=header)
            print(tabulate(df, headers='keys', tablefmt='psql'))
        else:
            raise NotImplementedError(format)

    def _process(self,segm:AnySegmentType) -> bool:
        """ """
        handled=False
        #print(segm.target,segm.start,segm.end,segm.length)
        desires_new_position = False
        if self.min_nt_length and segm.length < self.min_nt_length:
            # explicit request to fully skip this tiny segment
            # TODO: place is some kind of queue to report on elsewhere?
            return False
        elif self.current_chrom == None:
            self._place_segment_on_new_row_hook(segm)
            handled=True
        elif self.current_chrom != segm.target:
            self._flush_current_chrom_phase_table()
            self.last_pos = 0
            self._place_segment_on_new_row_hook(segm)
            handled=True
        elif self.last_pos <= segm.start:
            self._place_segment_on_new_row_hook(segm)
            handled=True
        elif self.last_pos >= segm.end:
            # included
            desires_new_position = False
        elif self.max_overlap_rel and self.current_segment._get_relative_overlap(segm) > self.max_overlap_rel:
            desires_new_position = False
        elif self.max_overlap_abs and self.current_segment._get_absolute_overlap(segm) > self.max_overlap_abs:
            desires_new_position = False
        elif self.min_unique_contribution and segm.end - self.last_pos >= self.min_unique_contribution:
            desires_new_position = True
        elif self.min_unique_contribution and (segm.start - self.current_segment.start) >= self.min_unique_contribution:
            desires_new_position = True
        else:
            desires_new_position = True

        if not handled:
            if desires_new_position:
                self._place_segment_on_new_row_hook(segm)
                handled = True
            else:
                self._place_segment_on_existing_row_hook(segm)
                handled = True

        # eof _process_segment
        self._post_process_segment(segm)
        return handled


class PySamInputDataProcessor(InputDataProcessor):
    """ """
    def _expose_inputdata(self,samfile) -> Any:
        """ """
        if isinstance(samfile,pysam.AlignmentFile):
            sam = samfile
        elif isinstance(samfile,Path):
            sam = pysam.AlignmentFile(samfile, "rb")
        elif type(samfile) == str and os.path.isfile(samfile):
            sam = pysam.AlignmentFile(samfile, "rb")
        else:
            raise NotImplementedError(samfile)
        return sam

    def _loop_over_inputdata(self):
        """ """
        with self.inputdata:
            mapping_idx = 0
            for read in self.inputdata:
                mapping_idx += 1
                if read.is_unmapped:
                    continue
                self._process(self._wrap(read, mapping_idx))

    def _wrap(self, read:pysam.AlignedSegment, idx:int) -> AnySegmentType:
        """ """
        ref = GenericSegment(
            read.reference_name, read.reference_start, read.reference_end,read.query_name,
            "-" if read.is_reverse else "+", 0
        )
        qry = GenericSegment(
            read.query_name, read.query_alignment_start, read.query_alignment_end,read.reference_name,
            "+", read.qlen)
        return AlignedSegment(
            ref=ref,
            query=qry,
            score=float(read.mapping_quality),
            idx=idx,
            raw=read,
        )


class FilteredBamOutputDataProcessor(OutputDataProcessor):

    def _prepare_output_files(self):
        """ hook which opens the provided output BAM file """
        outfname = self.output_config['outbam']
        self.outfile = pysam.AlignmentFile(outfname, "wb", template=self.inputdata)

    def _flush_scaffolded_input_data(self, out=Any):
        """ flush the SAM seqments from the primary phases, ditch the "bumped" segments """
        # we want to maintain input ORDER (segm.idx), to maintain possibly ordered BAM
        original_ordering_data = []
        for rowId in range(len(self.active_phases[0])):
            if self.active_phases[0][rowId].target != self.current_chrom:
                continue
            # loop over the active non-bumped phases and add "X1:Z:PhaseId_\d" tag
            for phaseId in range(1,self.ploidy+1):
                segm = self.active_phases[(phaseId-1)][rowId]
                if not segm:
                    continue
                # get (raw input) read and manipulate its tags
                read = segm.raw
                tags = read.get_tags()
                tags.append(("X1","PhaseId_%d"%phaseId))
                read.set_tags(tags)
                # store order information (maintains sorted BAM)
                row = ((segm.idx,rowId,(phaseId-1)))
                original_ordering_data.append(row)

        # now, we have to maintain ORDER (segm.idx), to maintain original order
        original_ordering_data.sort()
        for idx,rowId,phaseIdZb in original_ordering_data:
            read = self.active_phases[phaseIdZb][rowId].raw
            self.outfile.write(read)

class WgaScaffolderFromPysam(MatrixScaffolder,PySamInputDataProcessor):
    """ InputDataProcessor switched to SAM/BAM input from file using pysam """
    pass

class WgaScaffolderFilteredPysam(WgaScaffolderFromPysam,FilteredBamOutputDataProcessor):
    """ OutputDataProcessor switched to writing SAM/BAM, but excluding bumped segments (using pysam) """
    pass


if __name__ == "__main__":
    MSD = MatrixScaffolderDefaults
    sys.stderr.write("#"+str(MSD)+"\n")
    MSD.max_overlap_abs = 500
    MSD.min_unique_contribution = 1000

    try:
        scenario = int(sys.argv[2])
    except (IndexError, ValueError) as e:
        scenario = 1

    if scenario == 2:
        MSD.ploidy = 2
        MSD.max_overlap_abs = 10000
        #MSD.max_overlap_abs = 500
        MSD.min_unique_contribution = 25000

        MSD.min_nt_length = 5000

    elif scenario == 3:
        #MSD.max_overlap_abs = 1000
        MSD.max_overlap_abs = 5000
        MSD.min_unique_contribution = 200
        MSD.min_nt_length = 1500

    elif scenario == 4:
        MSD.ploidy = 4
        MSD.max_overlap_abs = 10000
        #MSD.max_overlap_abs = 500
        MSD.min_unique_contribution = 25000
        MSD.min_nt_length = 7500



    sys.stderr.write("#"+str(MSD)+"\n")
    scf = WgaScaffolderFromPysam(**MSD.__dict__)
    scf.run(sys.argv[1])
    if '--tsv' in sys.argv:
        scf.report(format='tsv')
    else:
        scf.report()

    ## with configured output
    scf = WgaScaffolderFilteredPysam(**MSD.__dict__).configure_output(outbam="/tmp/output.bam")
    scf.run(sys.argv[1])
    if '--tsv' in sys.argv:
        scf.report(format='tsv')
    else:
        scf.report()

"""
samtools view -h --remove-flags 0x900 /tmp/output.bam | ~/software/k8-1.2/k8-x86_64-Linux /home/avdb/software/minimap2/misc/paftools.js  sam2paf - | sort -k6,6 -k8,8n > /tmp/output.sorted.paf
/home/avdb/software/k8-1.2/k8-x86_64-Linux /home/avdb/software/minimap2/misc/paftools.js call -f data/cladegroup-cauris/GCA_002759435.3__external__CM105539.1.fa /tmp/output.sorted.paf > cladeIV-on-CWZ.vcf 

grep -v "^#" cladeIV.vcf | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $0}' | cut -f 1-11 | bedtools cluster -d 15
"""

