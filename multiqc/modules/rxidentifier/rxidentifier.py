import logging
from collections import OrderedDict, defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph #TODO is this what we want?

# Initialize the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """RxIdentifier module"""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="RxIdentifier",
            anchor="rxidentifier",
            href="https://github.com/genome/qc-analysis-pipeline",
            info="implementation of Rx algorithm in a WDL workflow",
            doi="https://doi.org/10.1371/journal.pone.0163019.s003",
        )

        # Find and load any RxIdentifier outputs
        self.rxidentifier_data = dict()
        for f in self.find_log_files("rxidentifier"):
            self.parse_rxidentifier(f)

        # Filter to strip out ignored sample names
        self.rxidentifer_data = self.ignore_samples(self.rxidentifier_data)

        # Warning when no files are found
        if len(self.rxidentifier_data) == 0:
            raise UserWarning

        # Write parsed data to a file
        self.write_data_file(self.rxidentifier_data, "multiqc_rxidentiifer")

        # Basic Stats Table
        self.rxidentifier_general_stats_table()

        # Rx Plot
        self.rxidentifier_plot()

    def parse_rxidentifier(self, f):
        line = f["f"]
        (sample_id, total_map, rx_score, confinterval_size, ci_text, result) = line.split("\t")
        result = result.strip()
        sample_id = self.clean_s_name(sample_id, f)

        if sample_id in self.rxidentifier_data:
            log.debug(f"Duplicate sample name found! Overwriting: {sample_id}")

        self.add_data_source(f, sample_id)
        self.rxidentifier_data[sample_id] = {
            "total_map": int(total_map),
            "rx_score": float(rx_score),
            "confinterval_size": confinterval_size,
            "ci_text": ci_text,
            "ci_lower": float(rx_score)-float(confinterval_size),
            "ci_upper": float(rx_score)+float(confinterval_size),
            "rx_result": result,
        }

    def rxidentifier_general_stats_table(self):
        """Take the parsed stats from the RxIdentifier results and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["total_map"] = {
            "title": "Total Mapped Reads",
            "description": "Total Mapped Reads",
            "min": 0,
            "format": '{:,.0f}',
            "scale": "Greens",
        }
        headers["rx_score"] = {
            "title": "Rx Score",
            "description": "Numerical Rx Value",
            "min": 0,
            "max": 1,
            "format": '{:,.3f}',
            "scale": "RdBu",
        }
        headers["rx_result"] = {
            "title": "Rx Result",
            "description": "Evaluated Rx Value",
            "scale": "Set1",
        }

        general_data = dict()
        for sample in self.rxidentifier_data.keys():
            general_data[sample] = {
                "total_map": self.rxidentifier_data[sample]["total_map"],
                "rx_score": self.rxidentifier_data[sample]["rx_score"],
                "rx_result": self.rxidentifier_data[sample]["rx_result"],
            }

        self.general_stats_addcols(general_data, headers)

    def rxidentifier_plot(self):
        #
        cats = OrderedDict()
        cats["rx_score"] = {
            "name": 'Rx score',
        }

        scoreplot_config = {
            "id": "rxidentifier_rx_score_plot",
            "title": "RxIdentifier: Rx score by sample",
            "ylab": "Rx",
            "ymax": 1.0,
            "tt_label": "{x}: {y:.8f}",
            "cpswitch": False,
            "tt_decimals": 8,
        }

#        self.add_section(
#            name="Rx values",
#            anchor="rxidentifier_rx_values",
#            description="""
#            This plot shows the reported Rx values for each sample
#            """,
#            plot=bargraph.plot(self.rxidentifier_data, cats, pconfig=scoreplot_config)
#        )

        cats = OrderedDict()
        cats["XX"] = {
            "name": "XX",
            "color": "#2b83ba",
        }
        cats["consistent with XX but not XY"] = {
            "name": "consistent with XX but not XY",
            "color": "#abdda4",
        }
        cats["Not Assigned"] = {
            "name": "Not Assigned",
            "color": "#ffffbf",
        }
        cats["consistent with XY but not XX"] = {
            "name": "consistent with XY but not XX",
            "color": "#fdae61",
        }
        cats["XY"] = {
            "name": "XY",
            "color": "#d7191c",
        }

        scoreplot_data = {}
        for sample in self.rxidentifer_data.keys():
            sample_score = self.rxidentifier_data[sample]["rx_score"]
            sample_cat = self.rxidentifer_data[sample]["rx_result"]

            scoreplot_data[sample] = {sample_cat: sample_score}

        self.add_section(
            name="Rx values",
            anchor="rxidentifier_rx_values",
            description="""
            This plot shows the reported Rx values for each sample
            """,
            plot=bargraph.plot(scoreplot_data, cats, pconfig=scoreplot_config)
        )

        classifications = defaultdict(int)
        for data in self.rxidentifier_data.values():
            classifications[data["rx_result"]] += 1

        resultplot_config = {
            "id": "rxidentifier_rx_result_plot",
            "title": "RxIdentifier: Samples per report result",
        }

        self.add_section(
            name="Rx Results",
            anchor="rxidentifier_rx_results",
            description="""
            This plot shows the number of samples reported for each result
            """,
            plot=bargraph.plot({"Rx results": classifications}, cats, pconfig=resultplot_config)
        )

