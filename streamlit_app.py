import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import boto3
# from dotenv import load_dotenv

# Load environment variables from .env file if needed
# load_dotenv()
# st.balloons()

# Read AWS credentials from environment variables
# @st.cache_data
# def get_adata():
#     aws_access_key_id = st.secrets['AWS_ACCESS_KEY_ID']
#     aws_secret_access_key = st.secrets['AWS_SECRET_ACCESS_KEY']

#     # Initialize S3 client
#     s3 = boto3.client('s3', aws_access_key_id=aws_access_key_id, aws_secret_access_key=aws_secret_access_key)

#     bucket_name = 'testzl57208'
#     file_key = 'poplar_axl_bud_A_adata.h5ad'
#     local_file_name = 'poplar_axl_bud_A_adata.h5ad'
#     # print(f'Reading file: {local_file_name}')

#     # Download file from S3
#     s3.download_file(bucket_name, file_key, local_file_name)

#     # Read the AnnData file using Scanpy
#     adata1 = sc.read(local_file_name)
#     return adata1



def main():
    ## dic for matching input sample same and stored adata file name
    sample_dic = {
        'scRNA_Leaf_1': 'sobj_Leaf_1_slim.h5ad',
        'scRNA_Leaf_3': 'sobj_Leaf_3_slim.h5ad',
        'scRNA_Leaf_5': 'sobj_Leaf_5_slim.h5ad',
        'scRNA_Primary_stem': 'sobj_Primary_stem_slim.h5ad',
        'scRNA_Secondary_stem': 'sobj_Secondary_stem_slim.h5ad',
        'spRNA_SAM': 'sobj_SAM_slim.h5ad'
    }
    
    # Organize tissues by library type
    tissues_by_lib = {
        'scRNA': ['Leaf_1', 'Leaf_3', 'Leaf_5', 'Primary_stem', 'Secondary_stem'],
        'spRNA': ['SAM']
    }
        
    ###############################################
    ##                 Sidebar                  ###
    ###############################################
    # use from() if want to add a Submit button before change everything
    # with st.sidebar.form():           
    
    st.sidebar.image('images/BioPoplar_Logo2.png', width="stretch")
    st.sidebar.header('Plot Configuration')
    st.sidebar.markdown('## Please select a dataset:')
    lib_type = st.sidebar.selectbox('Library Type', ['---Please choose---','scRNA'])    # 'spRNA'
    
    if lib_type == 'scRNA':
        tissue = st.sidebar.selectbox('Tissue', ['---Please choose---'] + tissues_by_lib['scRNA'])
    elif lib_type == 'spRNA':
        tissue = st.sidebar.selectbox('Tissue', ['---Please choose---'] + tissues_by_lib['spRNA'])
    else:
        tissue = None
    
    adata = None
    ## sameple name based on lib_type and tissue selection
    
    sample_name = None
    gene_ids = []
    if lib_type != '---Please choose---' and tissue != '---Please choose---' and tissue != None:
        sample_name = lib_type + '_' + tissue
        # st.sidebar.markdown('**Selected data type:**  \n%s' % sample_name)
        if sample_name not in sample_dic.keys():
            sample_name = None
            st.sidebar.write('No data available. Please select another dataset.')
        else:
            filename = sample_dic[sample_name]
            # st.sidebar.markdown('**Selected dataset:**  \n%s' % filename)  
            adatafile = 'data/' + filename
            # st.sidebar.markdown('**adata file:**  \n%s' % adatafile)  
            adata = sc.read_h5ad(adatafile)
        
        # ## for testing... remove later
        # adata = sc.read_h5ad('data/scrna_poplar_shoot_merged_0.9.h5ad')
        
        # read the genes
        if adata is not None:
            gene_ids = adata.var.index.tolist()
    else:
        st.sidebar.write('Please select a data to explore the genes')
    # print(filename)

    # gene_ids = ['test','ann1.Glyma.02G228100', 'ann1.Glyma.15G127900']
    
    plot_type = None
    if sample_name != None and sample_name != '---Please choose---':
        st.sidebar.markdown('## Please select gene to plot:')
        if lib_type == 'scRNA':
            plot_type = st.sidebar.radio('Select plot type', ['UMAP', 'Dotplot', 'Heatmap'], horizontal=True)
        elif lib_type == 'spRNA':
            plot_type = st.sidebar.radio('Select plot type', ['Spatial UMAP', 'Spatial', 'Dotplot', 'Heatmap'], horizontal=True)
        
        # Select gene(s) from available gene IDs using a multiselect
        if gene_ids:
            gene_list = st.sidebar.multiselect('Select gene(s)', options=gene_ids, default=[gene_ids[0]])
        else:
            gene_list = []
            st.sidebar.write('No genes available for selected dataset.')

        gene_name = gene_list[0] if gene_list else ''
        
            
    ###############################################
    ##                Main page                 ###
    ###############################################
    
    ## check if sample was selected:
    if sample_name and sample_name != '---Please choose---':

        
        # selected gene for plot
        variables_to_plot = ['seurat_clusters']        ## cell typs not for Violin plot
        if gene_name:
            variables_to_plot.append(gene_name)
        # st.markdown('Selected gene `%s`' % gene_name)
        
        if plot_type == 'UMAP':
            st.subheader('UMAP Plot')
            fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
            if len(variables_to_plot) == 1:
                axs = [axs]
            for ax, gene in zip(axs, variables_to_plot):
                sc.pl.umap(adata, color=gene, ax=ax, show=False)
                plt.subplots_adjust(wspace=1.2)
            st.pyplot(fig)
            plt.close(fig)
        elif plot_type == 'Spatial UMAP':
            st.subheader('Spatial UMAP Plot')
            fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
            if len(variables_to_plot) == 1:
                axs = [axs]
            for ax, gene in zip(axs, variables_to_plot):
                sc.pl.umap(adata, color=gene, ax=ax, show=False)
                plt.subplots_adjust(wspace=1.2)
            st.pyplot(fig)
            plt.close(fig)
        elif plot_type == 'Spatial':
            st.subheader('Spatial Plot: ' + sample_name)
            if lib_type != 'spRNA':
                st.error(' :crying_cat_face: Spatial Map not available. \nSelected data is not spatial transcriptomics')
            else:
                fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
                if len(variables_to_plot) == 1:
                    axs = [axs]
                for ax, gene in zip(axs, variables_to_plot):
                    sc.pl.spatial(adata, color=gene, ax=ax, show=False)
                    plt.subplots_adjust(wspace=1.2)
                st.pyplot(fig)
        elif plot_type == 'Dotplot':
            st.markdown('**Dotplot**')
            if gene_list:
                fig, ax = plt.subplots(figsize=(12, 6))
                sc.pl.dotplot(adata, gene_list, groupby='seurat_clusters', ax=ax)
                st.pyplot(fig)
                plt.close(fig)
            else:
                st.error(':point_left: Please enter gene names in the form')
        elif plot_type == 'Heatmap':
            st.markdown('**Heatmap**')
            if gene_list:
                hm = sc.pl.heatmap(
                adata, 
                gene_list, 
                groupby='seurat_clusters', 
                show=False, 
                swap_axes=True)
                fig = plt.gcf()
                fig.set_size_inches(12, 8)
                st.pyplot(fig)
                plt.close(fig)
            else:
                st.error(':point_left: Please enter gene names in the form')
                
    if plot_type == None:
        st.markdown('## Poplar single-cell Gene Viewer')
        st.markdown('Welcome to the poplar single-cell database. Please select the dataset and genes from the sidebar to start.')
        # st.markdown('[paper title]([paper link]), Giabardo et al., ')
        st.image('./images/scRNA_workflow.png')
    ## footnote
    # footnote = """
    #             <hr>
    #             <p style='font-size: small;'>For issues and questions, please contact 
    #             <a href='mailto:luoziliang@uga.edu'>Ziliang</a>.</p>
    #             """
    # st.markdown(footnote, unsafe_allow_html=True)
    
    
if __name__ == '__main__':
    main()