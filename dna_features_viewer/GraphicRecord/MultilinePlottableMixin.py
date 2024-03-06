import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy


class MultilinePlottableMixin:
    def plot_on_multiple_lines(
        self,
        n_lines=None,
        nucl_per_line=None,
        plot_sequence=False,
        figure_width="auto",
        translation_params=None,
        **plot_params
    ):
        """Plot the features on different lines (one Matplotlib ax per line).

        Parameters
        ----------

        n_lines
          Number of lines on which the record will be plotted. A number of
          nucleotides per line can be provided instead (see below).

        nucl_per_line
          Number of nucleotides to be represented on every line (determines
          the number of lines ``n_lines``).

        plot_sequence
          Whether to plot the nucleotide sequence on each line.

        figure_width
          Width of the figure in inches. Leave to auto for a width of either 10
          (if not sequence is plotted) or 0.15*nucl_per_line inches
          (if a sequence is plotted).

        translation_params
          Parameters for sequence translation. By default (``None``), it does
          not plot a translated sequence.

        **plot_params
          Parameters from ``graphic_record.plot()`` to be used in the plotting
          of the individual lines. This includes ``draw_line``, ``with_ruler``,
          ``annotate_inline``, ``plot_sequence``,
          ``evelate_outline_annotations``, ``strand_in_label_pixel_threshold``.

        Returns
        -------

        figure, axes
          The matplotlib figure and axes generated.
        """

        if n_lines is None:
            n_lines = int(numpy.ceil(self.sequence_length / nucl_per_line))
        else:
            nucl_per_line = self.sequence_length // n_lines + 1

        if figure_width == "auto":
            if plot_sequence:
                figure_width = 0.15 * nucl_per_line
            else:
                figure_width = 10

        figures_heights = []

        def plot_line(line_index, ax=None, translation_params=None):
            first, last = self.first_index, self.last_index
            line_start = first + line_index * nucl_per_line
            line_virtual_end = first + (line_index + 1) * nucl_per_line
            line_end = min(last, line_virtual_end)
            line_record = self.crop((line_start, line_end))
            line_ax, _ = line_record.plot(
                figure_width=figure_width,
                x_lim=(line_start, line_virtual_end),
                ax=ax,
                plot_sequence=plot_sequence,
                **plot_params
            )

            if translation_params is not None:
                line_record.plot_translation(ax=line_ax, **translation_params)

            return line_ax

        for line_index in range(n_lines):
            line_ax = plot_line(line_index, translation_params=translation_params)
            figures_heights.append(line_ax.figure.get_figheight())
            plt.close(line_ax.figure)
        fig, axes = plt.subplots(
            n_lines,
            1,
            gridspec_kw={"height_ratios": figures_heights},
            figsize=(figure_width, 0.9 * sum(figures_heights)),
        )
        if n_lines == 1:
            axes = [axes]
        for line_index, ax in enumerate(axes):
            plot_line(line_index, ax=ax, translation_params=translation_params)
        fig.tight_layout()
        return fig, axes
    
    def plot_on_multiple_lines_with_density(
        self,
        density_list : list,
        density_label : str = '',
        n_lines=None,
        nucl_per_line=None,
        plot_sequence=False,
        figure_width="auto",
        translation_params=None,
        **plot_params
    ):
        """Plot the features on different lines (one Matplotlib ax per line).
           with a density plot below each line

        Parameters
        ----------
        density_list
          list of density values

        n_lines
          Number of lines on which the record will be plotted. A number of
          nucleotides per line can be provided instead (see below).

        nucl_per_line
          Number of nucleotides to be represented on every line (determines
          the number of lines ``n_lines``).

        plot_sequence
          Whether to plot the nucleotide sequence on each line.

        figure_width
          Width of the figure in inches. Leave to auto for a width of either 10
          (if not sequence is plotted) or 0.15*nucl_per_line inches
          (if a sequence is plotted).

        translation_params
          Parameters for sequence translation. By default (``None``), it does
          not plot a translated sequence.

        **plot_params
          Parameters from ``graphic_record.plot()`` to be used in the plotting
          of the individual lines. This includes ``draw_line``, ``with_ruler``,
          ``annotate_inline``, ``plot_sequence``,
          ``evelate_outline_annotations``, ``strand_in_label_pixel_threshold``.

        Returns
        -------

        figure, axes
          The matplotlib figure and axes generated.
        """

        if n_lines is None:
            n_lines = int(numpy.ceil(self.sequence_length / nucl_per_line))
        else:
            nucl_per_line = self.sequence_length // n_lines + 1

        if figure_width == "auto":
            if plot_sequence:
                figure_width = 0.15 * nucl_per_line
            else:
                figure_width = 10

        figures_heights = []

        def plot_line(line_index, ax=None, translation_params=None):
            first, last = self.first_index, self.last_index
            line_start = first + line_index * nucl_per_line
            line_virtual_end = first + (line_index + 1) * nucl_per_line
            line_end = min(last, line_virtual_end)
            line_record = self.crop((line_start, line_end))
            line_ax, _ = line_record.plot(
                figure_width=figure_width,
                x_lim=(line_start, line_virtual_end),
                ax=ax,
                plot_sequence=plot_sequence,
                **plot_params
            )
            if translation_params is not None:
                line_record.plot_translation(ax=line_ax, **translation_params)

            return line_ax
        
        def plot_density(line_index,ax):
            """
            
            """
            first, last = self.first_index, self.last_index
            line_start = first + line_index * nucl_per_line
            line_virtual_end = first + (line_index + 1) * nucl_per_line
            # line_end = min(last, line_virtual_end)
            # xx = list(range(line_start,line_end))
            idx_start =  line_index * nucl_per_line
            idx_end = min(idx_start + nucl_per_line,len(density_list))
            yy = density_list[idx_start:idx_end]
            xx = numpy.arange(len(yy)) + line_start
            # print('xx',xx,'yy',yy)
            ax.fill_between(xx, yy, alpha=0.3)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(density_label)
            return ax

  
        for line_index in range(n_lines):
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(figure_width, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
            plot_line(line_index, translation_params=translation_params, ax=ax1)
            plot_density(line_index, ax2)
            figures_heights.append(fig.get_figheight())
            fig.show()
            plt.close(fig)

        fig = plt.figure(layout='constrained',
                          figsize=(figure_width, 0.9 * sum(figures_heights)))
        subfigs = fig.subfigures(n_lines, 1,
                                 height_ratios=figures_heights,
                                 width_ratios=[figure_width],
                                 hspace=0
                                 ) 

        if n_lines == 1:
          subfigs = [subfigs]

        for line_index, subfig in enumerate(subfigs):
            ax1, ax2 = subfig.subplots(2, 1, 
                                       sharex=True, 
                                       gridspec_kw={"height_ratios": [4, 1]}
                                       )
            plot_line(line_index, ax=ax1 ,translation_params=translation_params)
            plot_density(line_index, ax2)
    

        return fig

    def plot_on_multiple_pages(
        self,
        pdf_target,
        n_lines=None,
        nucl_per_line=None,
        lines_per_page=5,
        figure_width="auto",
        translation_params=None,
        **plot_params
    ):
        """Plot the features on different lines on different pages of a PDF.

        This function returns None

        Parameters
        ----------

        pdf_target
          Either a path to a PDF, or a file(-like) handle.

        n_lines
          Number of lines on which the record will be plotted. A number of
          nucleotides per line can be provided instead (see below).

        nucl_per_line
          Number of nucleotides to be represented on every line (determines
          the number of lines ``n_lines``).

        lines_per_page
          Number of lines on each page.

        plot_sequence
          Whether to plot the nucleotide sequence on each line.

        figure_width
          Width of the figure in inches. Leave to auto for a width of either 10
          (if not sequence is plotted) or 0.15*nucl_per_line inches
          (if a sequence is plotted).

        translation_params
          Parameters for sequence translation. By default (``None``), it does
          not plot a translated sequence.

        **plot_params
          Parameters from ``graphic_record.plot()`` to be used in the plotting
          of the individual lines. This includes ``draw_line``, ``with_ruler``,
          ``annotate_inline``, ``plot_sequence``,
          ``evelate_outline_annotations``, ``strand_in_label_pixel_threshold``.
        """
        nucl_per_page = nucl_per_line * lines_per_page
        number_of_pages = int(numpy.ceil(self.sequence_length / nucl_per_page))
        with PdfPages(pdf_target) as pdf:
            for page_index in range(number_of_pages):
                first, last = self.first_index, self.last_index
                page_start = first + page_index * nucl_per_page
                page_end = first + (page_index + 1) * nucl_per_page
                page_end = min(last, page_end)
                page_record = self.crop((page_start, page_end))
                fig, axes = page_record.plot_on_multiple_lines(
                    nucl_per_line=nucl_per_line,
                    figure_width=figure_width,
                    translation_params=translation_params,
                    **plot_params
                )
                pdf.savefig(fig)
                plt.close(fig)
