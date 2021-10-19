# Get the total excess precipitation before and after storm water reduction
import json
import pandas as pd
import pathlib as pl


def extract_event_totals(forcing_files: list) -> dict:
    """Given the paths to the forcing or unreduced forcing json, extract the total excess precipitation for each domain and event."""
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


outputs_dir = pl.Path("/Outputs")
forcing_dir = outputs_dir / "CedarRapids_P01_Forcing"
metadata_dir = outputs_dir / "CedarRapids_P01_Metadata"
event_excess_path = outputs_dir / "CedarRapids_Event_Excess_Precipitation.xlsx"
sheet_name_prefix = "P01"

forcing_files = list(forcing_dir.glob("**/*.json"))
unreduced_forcing_files = [file for file in metadata_dir.glob("**/*.json") if "Unreduced" in file.stem]
print(f"{len(forcing_files)} forcing files and {len(unreduced_forcing_files)} unreduced forcing files identified")

event_totals = extract_event_totals(forcing_files)
unreduced_event_totals = extract_event_totals(unreduced_forcing_files)

writer = pd.ExcelWriter(event_excess_path)
for domain in event_totals.keys():
    event_totals_df = create_event_total_dataframe(event_totals, unreduced_event_totals, domain)
    event_totals_df.to_excel(writer, sheet_name=f"{sheet_name_prefix}_{domain}")
writer.save()
