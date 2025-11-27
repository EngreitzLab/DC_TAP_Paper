import pandas as pd


def merge_singleton(df_s3, df_sing):
    """
    Helper method to merge singleton dataframe to table_s3 dataframe
    """
    # Split comma-separated guide_ids into lists
    df_s3["guide_list"] = df_s3["guide_ids"].str.split(",")

    # Explode to one row per guid
    df_s3_expanded = df_s3.explode("guide_list")
    df_s3_expanded["guide_list"] = df_s3_expanded["guide_list"].str.strip()

    # Merge with df_sing to get effect_size for each guide
    merged = df_s3_expanded.merge(
        df_sing[
            [
                "grna_id",
                "response_id",
                "pct_change_effect_size",
                "n_nonzero_trt",
                "n_nonzero_cntrl",
                "standard_error_pct_change",
            ]
        ],
        left_on=["guide_list", "gene_id"],
        right_on=["grna_id", "response_id"],
        how="left",
    )

    # Re-label column names
    merged = merged.rename(
        columns={
            "response_id": "singleton_response_id",
            "pct_change_effect_size_x": "pct_change_effect_size",  # fix original colname
            "pct_change_effect_size_y": "singleton_pct_change_effect_size",
            "n_nonzero_trt": "singleton_n_nonzero_trt",
            "n_nonzero_cntrl": "singleton_n_nonzero_cntrl",
            "standard_error_pct_change_x": "standard_error_pct_change",  # fix original colname
            "standard_error_pct_change_y": "singleton_standard_error_pct_change",
        }
    )

    return merged
