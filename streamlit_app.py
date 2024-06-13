import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3
# from dotenv import load_dotenv

# Load environment variables from .env file if needed
# load_dotenv()
st.balloons()

# Read AWS credentials from environment variables
aws_access_key_id = st.secrets['AWS_ACCESS_KEY_ID']
aws_secret_access_key = st.secrets['AWS_SECRET_ACCESS_KEY']

# Initialize S3 client
s3 = boto3.client('s3', aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key)

bucket_name = 'testzl57208'
file_key = 'poplar_axl_bud_A_adata.h5ad'
local_file_name = 'poplar_axl_bud_A_adata.h5ad'

# Download file from S3
s3.download_file(bucket_name, file_key, local_file_name)

# Read the AnnData file using Scanpy
adata1 = sc.read(local_file_name)
# adata1 = sc.read_h5ad('')

# Get the list of gene IDs
gene_ids = adata1.var.index.tolist()

# Sidebar for plot type selection and gene name input
st.sidebar.header('Plot Configuration')
st.write('Please select the plot type (UMAP or spatial), and gene ID for the plot')
plot_type = st.sidebar.selectbox('Select plot type', ['UMAP', 'Spatial'])
gene_name = st.sidebar.selectbox('Enter gene name for expression plot', ['', *gene_ids])

# Main panel
st.title('Spatial Transcriptome Visualization')

variables_to_plot = ['clusters']
if gene_name:
    variables_to_plot.append(gene_name)

if plot_type == 'UMAP':
    st.subheader('UMAP Plot')
    fig, axs = plt.subplots(1, len(variables_to_plot), figsize=(5 * len(variables_to_plot), 5))
    if len(variables_to_plot) == 1:
        axs = [axs]
    for ax, gene in zip(axs, variables_to_plot):
        sc.pl.umap(adata1, color=gene, ax=ax, show=False, wspace=0.6)
    st.pyplot(fig)
elif plot_type == 'Spatial':
    st.subheader('Spatial Plot')
    fig, axs = plt.subplots(1, len(variables_to_plot), figsize=(5 * len(variables_to_plot), 5))
    if len(variables_to_plot) == 1:
        axs = [axs]
    for ax, gene in zip(axs, variables_to_plot):
        sc.pl.spatial(adata1, color=gene, ax=ax, show=False, wspace=0.6)
    st.pyplot(fig)
else:
    st.write("Select a plot type from the sidebar to begin.")
