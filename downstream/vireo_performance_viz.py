        ##
        
        # UMAPs
        # labels_order = labels.value_counts().sort_values(ascending=False).index
        # gt_order = ground_truth.value_counts().sort_values(ascending=False).index
        # df_ = embs.join(a.obs)
        # df_['labels'] = pd.Categorical(labels, categories=labels_order)
        # df_['ground truth'] = pd.Categorical(ground_truth, categories=gt_order)
# 
        # # Show
        # fig, axs = plt.subplots(2,1, figsize=(5, 7), sharex=True)
        # draw_embeddings(
        #     df_, 
        #     cat='ground truth', 
        #     ax=axs[0], 
        #     legend_kwargs={'bbox_to_anchor' : (1,1), 'loc' : 'upper left'}
        # )
        # #axs[0].axis('off')
        # draw_embeddings(
        #     df_, 
        #     cat='labels', 
        #     title=f'vireoSNP: ARI {ARI:.2f}, NMI {NMI:.2f})',
        #     ax=axs[1], 
        #     legend_kwargs={'bbox_to_anchor' : (1,1), 'loc' : 'upper left'}
        # )
        # 
        # #axs[1].axis('off')
        # fig.tight_layout()
        # plt.show()

        # # Model fitting diagnostics 1
        # import matplotlib
        # matplotlib.use('macOSX')
        # 
        # fig = plt.figure(figsize=(11, 4))
        # plt.subplot(1, 2, 1)
        # plt.hist(_model.ELBO_inits)
        # plt.ylabel("Frequency")
        # plt.xlabel("ELBO in multiple initializations")
        # 
        # plt.subplot(1, 2, 2)
        # plt.plot(_model.ELBO_iters)
        # plt.xlabel("Iterations")
        # plt.ylabel("ELBO in a single initialization")
        # 
        # plt.tight_layout()
        # plt.show()

        # Visualize output
        # raw_col = cm.get_cmap('pink_r', 200)
        # new_col = np.vstack((raw_col(np.linspace(0, 0.7, 10)),
        #                      raw_col(np.linspace(0.7, 1, 90))))
        # segpink = ListedColormap(new_col, name='segpink')
# 
        # # Clonal assignment visualization
        # fig = plt.figure(figsize=(7, 4))
        # plt.subplot(1, 2, 1)
        # im = heat_matrix(_model.ID_prob, cmap="Blues", alpha=0.8,
        #                  display_value=False, row_sort=True)
        # plt.colorbar(im, fraction=0.046, pad=0.04)
        # plt.title("Assignment probability")
        # plt.xlabel("Clone")
        # plt.ylabel("%d cells" %(_model.n_cell))
        # plt.xticks(range(_model.n_donor))
# 
        # plt.subplot(1, 2, 2)
        # im = heat_matrix(_model.beta_mu, cmap=segpink, alpha=0.8,
        #                  display_value=False, row_sort=True)
        # plt.colorbar(im, fraction=0.046, pad=0.04)
        # plt.title("Mean allelic ratio")
        # plt.xlabel("Clone")
        # plt.ylabel("%d SNPs" %(_model.n_var))
        # plt.xticks(range(_model.n_donor))
# 
        # # Save
        # plt.tight_layout()
        # 
        # plt.show()
        # fig.savefig(os.path.join(path_results, sample, 'heatmaps.png'))

        ##