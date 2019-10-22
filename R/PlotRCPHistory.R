#' Plot the search path from a recursive clade partition search
#'
#' @description This function produces a list of tree-plots representing the
#' search-path during a recrusive clade partition search. Each plot has a label
#' on top showing:
#' \describe{
#'  \item{}{a number in parentheses (i) describes iteration i of the main loop.}
#'  \item{}{the score, the log-likelihood and the number of parameters of the model.}
#'  \item{}{a coloured node with a number i is the partition root for the iteration.}
#'  \item{}{Nodes in grey represent the potential shift points - these are
#'  descendants from the partition root, which have not been "cut out" by a shift
#'  and have at least q descendants, themselves.}
#'  \item{}{Letters in braces denote the candidate model-types for each shift-node.}
#' }
#' @param fit an object of S3 class `PCMFitModelMappings` returned by a call
#'   to `PCMFitMixed`.
#' @param sizeGreyNodepoints,sizeColorNodepoints,sizeBlackAllowedModelTypes,sizeColorAllowedModelTypes,sizeRankInQueue,vjustBlackAllowed,vjustColorAllowedModelTypes graphical parameters (see function description).
#' @param ... additional parameters passed to \code{\link{PCMTreePlot}}.
#' @return a list of annotated ggtree plots. Some of the entries in this list
#' can be NULL to indicate that no score improvement has been achieved at the
#' corresponding iteration. The example below shows how to filter these out.
#'
#' @importFrom data.table data.table setkey
#' @import PCMBase
#' @import PCMFit
#' @import ggtree
#' @import ggplot2
#'
#' @examples
#' lstPlots <- PlotRCPHistory(
#'   fitMappings_MGPM_B_best_clade_2_DataWithSEs, layout = "fan")
#' cowplot::plot_grid(plotlist = lstPlots[!sapply(lstPlots, is.null)])
#'
#' @export
PlotRCPHistory <- function(
  fit,
  sizeGreyNodepoints = 2.2, sizeColorNodepoints = 2.2,
  sizeBlackAllowedModelTypes = 1.4, sizeColorAllowedModelTypes = 1.4, sizeRankInQueue = 1.4,
  vjustBlackAllowedModelTypes = -1.6, vjustColorAllowedModelTypes = -1.6,
  ...) {
  tree <- PCMTree(fit$tree)

  treeRootInt <- PCMTreeNumTips(tree) + 1L
  PCMTreeSetLabels(tree)

  # need to create the color palette first.

  partNodes <- unique(unlist(lapply(seq_len(length(fit$mainLoopHistory)), function(i) {
    historyEntry <- fit$mainLoopHistory[[i]]
    historyEntry$headQPR_Partition
  })))

  colorsForPartNodes <- PCMColorPalette(length(partNodes), names = as.character(partNodes))

  plotList <- lapply(seq_len(length(fit$mainLoopHistory)), function(i) {
    #cat("Step ", i, "\n")

    historyEntry <- fit$mainLoopHistory[[i]]
    rootNodei <- fit$queuePartitionRoots[i, node]

    PCMTreeSetPartition(tree, historyEntry$headQPR_Partition)

    if(length(historyEntry$listPartitions) > 0) {
      dtCladePartition <- data.table(node=unique(unlist(historyEntry$listPartitions)))
      #print(dtCladePartition[, node])
      setkey(dtCladePartition, node)

      dtCladePartition[node%in%historyEntry$headQPR_Partition, selected:=TRUE]
      dtCladePartition[!(node%in%historyEntry$headQPR_Partition), candidate:=TRUE]

      # trying to reconstruct the remaining queue is hard - better to save it at runtime;
      remainingQueuei <- fit$queuePartitionRoots[i:historyEntry$lengthQPR, list(node = as.character(node))]

      remainingQueuei[, rankInQueue:=(i + .I - 1)]
      remainingQueuei <- remainingQueuei[!is.na(node)]
      setkey(remainingQueuei, node)

      dtCladePartition[, allowedModelTypes:=sapply(node, function(n) {
        iLabel <- as.integer(n)
        # we need the if(), because PCMTreeGetPartsForNodes returns an empty
        # vector for the root-node
        iRegime <- if(iLabel == treeRootInt) {
          historyEntry$headQPR_MappingIdx[1]
        } else {
          historyEntry$headQPR_MappingIdx[
            PCMTreeGetPartsForNodes(tree, iLabel)]
        }

        text <- do.call(
          paste,
          c(as.list(LETTERS[
            unique(c(iRegime, historyEntry$listAllowedModelTypesIndices[[n]]))]),
            list(sep="")))
        paste0("{",text,"}")
      })]

      dtCladePartition[
        node%in%historyEntry$headQPR_Partition,
        allowedModelTypesSelected:=allowedModelTypes]
      dtCladePartition[
        !(node%in%historyEntry$headQPR_Partition),
        allowedModelTypesCandidate:=allowedModelTypes]

      fitTablei <- RetrieveFittedModelsFromFitVectors(
        fit,
        LookupFit(
          tree = tree,
          modelTypes = fit$arguments$modelTypes,
          modelMapping = historyEntry$headQPR_Mapping,
          tableFits = fit$tableFits), setAttributes = TRUE)

      palette <- colorsForPartNodes[as.character(historyEntry$headQPR_Partition)]
      names(palette) <- as.character(PCMTreeGetRegimesForNodes(tree, nodes = as.integer(historyEntry$headQPR_Partition)))

      ploti <- PCMTreePlot(tree, palette = palette, ...) %<+% as.data.frame(dtCladePartition) %<+% as.data.frame(remainingQueuei) +
        geom_nodepoint(aes(shape=selected), size=sizeColorNodepoints, na.rm = TRUE) +
        geom_nodepoint(aes(shape=candidate), size=sizeGreyNodepoints, color = "grey", na.rm = TRUE) +
        geom_text(aes(label=allowedModelTypesSelected), size=sizeColorAllowedModelTypes, vjust=vjustColorAllowedModelTypes) +
        geom_text(aes(label=allowedModelTypesCandidate), color = "black", size=sizeBlackAllowedModelTypes, vjust=vjustBlackAllowedModelTypes) +
        geom_text(aes(label=rankInQueue), color="black", size=sizeRankInQueue) +
        ggtitle(paste0("(",i,") score=", round(fitTablei$score[[1]]), ", logLik=", round(fitTablei$logLik[[1]]), ", p=", fitTablei$df[[1]]))

      ploti
    } else {
      NULL
    }
  })
  plotList
}


