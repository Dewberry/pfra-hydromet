## Get the median precipitation and total excess precipitation before and after storm water reduction:
import json
import numpy as np
import pandas as pd
import pathlib as pl


def extract_median_precipitation(metadata_files: list) -> dict:
    """Given the paths to the metadata json, extract the median random precipitation for each grouped event."""
    event_totals = {}
    for metadata_file in metadata_files:
        domain = metadata_file.stem.split("_")[2]
        event_totals[domain] = {}
        with open(metadata_file) as f:
            metadata = json.load(f)
        for duration_metadata in metadata.values():
            groups = duration_metadata["BCName"][domain]["groups"]
            events_metadata = duration_metadata["BCName"][domain]["events_metadata"]
            events_metadata_df = pd.DataFrame.from_dict(events_metadata).set_index("EventID")
            for event, group in groups.items():
                event_totals[domain][event] = np.median(events_metadata_df.loc[group, "Random Precipitation"])
    return event_totals


def create_median_precipitation_dataframe(domain_precipitation: dict) -> pd.DataFrame:
    """Create a dataframe of the median random precipitation and grouped event ID."""
    median_precipitation_df = pd.DataFrame.from_dict({"Median_P": domain_precipitation})
    median_precipitation_df.sort_index(inplace=True)
    median_precipitation_df.index.name = "Event_ID"
    return median_precipitation_df


def extract_total_excess(forcing_files: list) -> dict:
    """Given the paths to the forcing or unreduced forcing json, extract the total excess precipitation for each domain
    and event."""
    event_totals = {}
    for forcing_file in forcing_files:
        with open(forcing_file) as f:
            event_data = json.load(f)
        domain = forcing_file.stem.split("_")[2]
        event_totals[domain] = {}
        for duration in list(event_data.keys()):
            for event, incremental_volume in event_data[duration]["BCName"][domain].items():
                event_totals[domain][event] = sum(incremental_volume)
    return event_totals


def create_event_total_dataframe(event_totals: dict, unreduced_event_totals: dict, domain: str) -> pd.DataFrame:
    """Create a dataframe of the excess precipitation and unreduced excess precipitation (if applicable) for the
    specified domain."""
    event_totals_df = pd.DataFrame.from_dict({"Excess_P": event_totals[domain]})
    if domain in unreduced_event_totals.keys():
        unreduced_event_totals_df = pd.DataFrame.from_dict({"Unreduced_Excess_P": unreduced_event_totals[domain]})
        if not all(unreduced_event_totals_df.index == event_totals_df.index):
            raise ValueError(
                "The excess precipitation and unreduced event IDs do not match, cannot concatenate the dataframes"
            )
        event_totals_df = pd.concat([unreduced_event_totals_df, event_totals_df], axis=1)
        if not all(event_totals_df["Excess_P"] <= event_totals_df["Unreduced_Excess_P"]):
            raise ValueError("One or more events have excess precipitation that is greater than the unreduced.")
    event_totals_df.sort_index(inplace=True)
    event_totals_df.index.name = "Event_ID"
    return event_totals_df


## Specify Paths and Variables:
outputs_dir = pl.Path("/Cedar_Rapids")
forcing_dir = outputs_dir / "CedarRapids_P01_Forcing"
metadata_dir = outputs_dir / "CedarRapids_P01_Metadata"
median_precip_path = outputs_dir / "CedarRapids_Median_Precipitation.xlsx"
event_excess_path = outputs_dir / "CedarRapids_Event_Excess_Precipitation.xlsx"
sheet_name_prefix = "P01"

## Identify Data Files:
forcing_files = list(forcing_dir.glob("**/*.json"))
unreduced_forcing_files = [file for file in metadata_dir.glob("**/*.json") if "Unreduced" in file.stem]
metadata_files = [file for file in metadata_dir.glob("**/*.json") if file.stem.split("_")[-1] == "Metadata"]
print(
    f"{len(forcing_files)} forcing files, {len(unreduced_forcing_files)} unreduced forcing, and {len(metadata_files)} files identified"
)


## Get the Median Precipitation:
median_precipitation = extract_median_precipitation(metadata_files)

writer = pd.ExcelWriter(median_precip_path)
for domain, domain_precipitation in median_precipitation.items():
    median_precipitation_df = create_median_precipitation_dataframe(domain_precipitation)
    median_precipitation_df.to_excel(writer, sheet_name=f"{sheet_name_prefix}_{domain}")
writer.save()


## Get the Total Excess Precipitation Before and After Storm Water Reduction
total_excess = extract_total_excess(forcing_files)
unreduced_total_excess = extract_total_excess(unreduced_forcing_files)

writer = pd.ExcelWriter(event_excess_path)
for domain in total_excess.keys():
    event_totals_df = create_event_total_dataframe(total_excess, unreduced_total_excess, domain)
    event_totals_df.to_excel(writer, sheet_name=f"{sheet_name_prefix}_{domain}")
writer.save()
