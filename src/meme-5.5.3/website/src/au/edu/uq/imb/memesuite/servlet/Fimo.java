package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.SequenceDB;
import au.edu.uq.imb.memesuite.db.SequencePrior;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

public class Fimo extends SubmitJob<Fimo.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentMotifs motifs;
  private ComponentSequences sequences;
  private ComponentBfile background;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions advancedOptions;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.fimo");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public MotifDataSource motifs;
    public SequenceInfo sequences;
    public SequencePrior priors;
    public Background background;
    public double pthresh;
    public boolean norc;
    public boolean pgc;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("motifs", motifs);
      out.property("sequences", sequences);
      out.property("priors", priors);
      out.property("background", background);
      out.property("pthresh", pthresh);
      out.property("norc", norc);
      out.property("pgc", pgc);
      out.endObject();
    }

    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;
    }

    @Override
    public boolean immediateRun() {
      return false;
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }
  
    @Override
    public String cmd() {
      StringBuilder args = new StringBuilder();
      if (sequences instanceof SequenceDataSource) {
        addArgs(args, "-upseqs", ((SequenceDataSource) sequences).getName());
      }
      if (sequences instanceof LociDataSource) {
        addArgs(args, "-bedfile");
        addArgs(args, "-genome", ((LociDataSource)sequences).getGenomeFileName());
      }
      if (priors != null) {
        addArgs(args, "-psp", priors.getPriorsName());
        addArgs(args, "-prior-dist", priors.getPriorsDistName());
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bfile", background.getBfile().getName());
      } else if (background.getSource() == Background.Source.NRDB) {
        addArgs(args, "-bfile", "--nrdb--"); 		// Use NRDB background
      } else if (background.getSource() == Background.Source.UNIFORM) {
        addArgs(args, "-bfile", "--uniform--");		// Use uniform background
      } else if (background.getSource() == Background.Source.MEME) {
        addArgs(args, "-bfile", "--motif--");		// Use background in motif file
      }
      addArgs(args, "-pvthresh", pthresh);
      if (norc) addArgs(args, "-norc");
      if (pgc) addArgs(args, "-parse-genomic-coord");
      addArgs(args, motifs.getName());
      if (sequences instanceof SequenceDB) {
        addArgs(args, ((SequenceDB) sequences).getSequenceName());
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (motifs != null) list.add(motifs);
      if (sequences != null && sequences instanceof SequenceDataSource) {
        list.add((SequenceDataSource) sequences);
      }
      if (background.getSource() == Background.Source.FILE) {
        list.add(background.getBfile());
      }
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (sequences != null && sequences instanceof SequenceDataSource) {
        if (!((SequenceDataSource) sequences).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((SequenceDataSource) sequences).getFile());
        }
      }
      if (motifs != null) {
        if (!motifs.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              motifs.getFile());
        }
      }
      if (background.getSource() == Background.Source.FILE) background.getBfile().getFile().delete();
    }
  }

  public Fimo() {
    super("FIMO", "FIMO");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/fimo.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/fimo_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
    jobDetails = new ComponentJobDetails(cache);
    advancedOptions = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("advanced_options"));
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
  }

  @Override
  public String title() {
    return tmplVerify.getSubtemplate("title").toString();
  }

  @Override
  public String subtitle() {
    return tmplVerify.getSubtemplate("subtitle").toString();
  }

  @Override
  public String logoPath() {
    return tmplVerify.getSubtemplate("logo").toString();
  }

  @Override
  public String logoAltText() {
    return tmplVerify.getSubtemplate("alt").toString();
  }

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), motifs.getHelp(),
        sequences.getHelp(), jobDetails.getHelp(), advancedOptions.getHelp(),
        submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("motifs", motifs.getComponent(request.getParameter("motifs_embed")));
    main.set("sequences", sequences.getComponent());
    main.set("bfile", background.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("advanced_options", advancedOptions.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback,
      HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    FileCoord.Name backgroundName = namer.createName("background");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the motifs
    data.motifs =  (MotifDataSource)motifs.getMotifs(motifsName, request, feedback);
    if (data.motifs != null) alph = data.motifs.getAlphabet();
    // get the sequences
    data.sequences = sequences.getSequences(alph, sequencesName, request, feedback);
    data.priors = sequences.getPrior(request);
    data.background = background.getBfile(backgroundName, request, feedback);
    // other options
    data.pthresh = WebUtils.paramNumber(feedback, "output p-value threshold",
      request, "output_pv", 0.0, 1.0, 1e-4);
    data.norc = WebUtils.paramBool(request, "norc");
    data.pgc = WebUtils.paramBool(request, "pgc");
    return data;
  }
}
