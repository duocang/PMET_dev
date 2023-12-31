package au.edu.uq.imb.memesuite.updatedb;

import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.util.MultiSourceStatus;
import au.edu.uq.imb.memesuite.util.Progress;
import static au.edu.uq.imb.memesuite.db.SQL.*;
import org.apache.commons.io.FileUtils;
import org.apache.commons.net.ftp.FTP;
import org.apache.commons.net.ftp.FTPClient;
import org.apache.commons.net.ftp.FTPReply;
import org.apache.commons.net.ftp.FTPFile;
import org.sqlite.SQLiteDataSource;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.SocketTimeoutException;
import java.net.URL;
import java.net.URLConnection;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.LinkedList;
import java.util.*;
import java.util.List;
import java.util.ListIterator;
import java.util.Properties;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantReadWriteLock;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
import org.json.simple.JSONValue;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.JSONParser;


/**
 * This is the base class of the sequence updaters. It implements shared functionality.
 */
public abstract class SequenceUpdater implements Callable<Void> {
  protected SQLiteDataSource dataSource;
  protected ReentrantReadWriteLock dbLock;
  protected File binDir;
  protected File dbDir;
  protected ExecutorService worker;
  protected CompletionService<Void> workerEcs;
  protected List<Future<Void>> workerTasks;
  protected Progress progress;
  protected MultiSourceStatus status;
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.updatedb.updater");
  private static final int CONNECT_TIMEOUT = (int)TimeUnit.SECONDS.toMillis(30);
  private static final int REPLY_TIMEOUT = (int)TimeUnit.SECONDS.toMillis(30);
  private static final long KEEP_ALIVE_TIMEOUT = TimeUnit.MINUTES.toMillis(4);

  protected static Pattern dnaFilePattern1 = Pattern.compile("^.*\\.dna\\.primary_assembly\\.fa\\.gz$");
  protected static Pattern dnaFilePattern2 = Pattern.compile("^.*\\.dna\\.toplevel\\.fa\\.gz$");
  protected static Pattern aaFilePattern = Pattern.compile("^.*\\.pep\\.all\\.fa\\.gz$");
  protected static Pattern aiFilePattern = Pattern.compile("^.*\\.pep\\.abinitio\\.fa\\.gz$");
  protected static Pattern releasePattern = Pattern.compile("/release-(\\d+)/");
  protected static Pattern collectionPattern = Pattern.compile("^.*collection$");
  protected static final int RETRY_COUNT = 3;
  protected static final int ERROR_COUNT = 10;

  public SequenceUpdater(String name, 
      SQLiteDataSource dataSource,
      ReentrantReadWriteLock dbLock, File binDir, File dbDir,
      ExecutorService worker, MultiSourceStatus status) {
    if (binDir == null) throw new NullPointerException("binDir is null");
    if (dbDir == null) throw new NullPointerException("dbDir is null");
    this.dataSource = dataSource;
    this.dbLock = dbLock;
    this.binDir = binDir;
    this.dbDir = dbDir;
    this.worker = worker;
    this.status = status;
    workerEcs = new ExecutorCompletionService<Void>(this.worker);
    workerTasks = new LinkedList<Future<Void>>();
    progress = new Progress(status, name);
  }

  // Read a JSON object from the Ensembl REST interface.
  protected String ensemblRestRead(URL url) {
    String output = "";
    Reader reader = null;

    try {
      URLConnection connection = url.openConnection();
      HttpURLConnection httpConnection = (HttpURLConnection)connection;

      httpConnection.setRequestProperty("Content-Type", "application/json");

      InputStream response = connection.getInputStream();
      int responseCode = httpConnection.getResponseCode();

      if (responseCode != 200) {
        throw new RuntimeException("Response code was not 200. Detected response was " + responseCode);
      }

      // Read from the REST interface.
      reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      StringBuilder builder = new StringBuilder();
      char[] buffer = new char[8192];
      int read;
      while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
        builder.append(buffer, 0, read);
      }
      output = builder.toString();
    } catch (IOException e) {
     logger.log(Level.SEVERE, "URL exception.");
    } finally {
      if (reader != null) try {
        reader.close();
      } catch (IOException logOrIgnore) {
        logOrIgnore.printStackTrace();
      }
    }

    return output;
  } // ensemblRestRead()

protected List<EnsemblGenome> ensemblAvailableGenomes(
  String genomeRelease, boolean abinitio, int retrieverId, String categoryName, 
  String host, String rootDir, String restSubdir, Map<String, String> name2Subdir) {
    progress.setTask("Querying available genomes using REST", 0, -1);
    List<EnsemblGenome> genomes = new ArrayList<EnsemblGenome>();

    // get a listing of species
    try {
      logger.log(Level.INFO, "Starting Update of " + categoryName);

      // Create the REST query.
      String server = "https://rest.ensembl.org";
      String ext = "/info/genomes/division/" + restSubdir + "?";
      URL url = new URL(server + ext);

      // Query the REST interface.
      String output = ensemblRestRead(url);

      // Parse the JSON string.
      JSONParser parser = new JSONParser();
      try {
        Object obj = parser.parse(output);
        JSONArray array = (JSONArray) obj;
        Iterator itr = array.iterator();
        boolean first_genome = true;
        while (itr.hasNext()) {
          JSONObject jenome = (JSONObject) itr.next();
          String displayName = jenome.get("display_name").toString();
          logger.log(Level.INFO, "Display Name: " + displayName);
          String scientificName = jenome.get("scientific_name").toString();
          String ftpDirName = jenome.get("name").toString();
          String path = rootDir;
          // Add the collection name to the path.
          String collection = name2Subdir.get(ftpDirName);
          if (collection != null) { path += collection; }
          path += ftpDirName;
          String dnaPath = path + "/dna/";
          String proteinPath = path + "/pep/";
          logger.log(Level.FINE, "Adding Ensembl DNA and protein genomes for " + displayName);
          if (! abinitio) {
            genomes.add(new EnsemblGenome(genomeRelease, retrieverId, categoryName, displayName, ftpDirName, AlphStd.DNA, false, dnaPath, host, dnaFilePattern1, dnaFilePattern2));
          }
          genomes.add(new EnsemblGenome(genomeRelease, retrieverId, categoryName, displayName, ftpDirName, AlphStd.PROTEIN, abinitio, proteinPath, host, abinitio ? aiFilePattern : aaFilePattern));
        }
      } catch(ParseException pe) {
        logger.log(Level.SEVERE, "position: " + pe.getPosition());
      }
    } catch (IOException e) {
     logger.log(Level.SEVERE, "URL exception in ensemblAvailableGenomes.");
    }
    return genomes;
  } // ensemblAvailableGenomes

  //
  // Make a Map of from the name of the species to its directory.
  // The species directory MUST be either:
  //   1) in the baseDir, or
  //   2) in a subdirectory of it named "*_collection".
  //
  protected Map<String, String> makeMapName2Subdir(
    FTPClient ftpClient, String rootDir
  ) throws IOException {
    Map<String, String> map = new HashMap<>();

    FTPFile[] topFiles = ftpClient.listFiles(rootDir);
    for (FTPFile aFile : topFiles) {
      Matcher m = collectionPattern.matcher(aFile.getName());
      if (m.find()) {
        String collectionName = aFile.getName();
        String subDir = rootDir + collectionName + "/";
        FTPFile[] subFiles = ftpClient.listFiles(subDir);
        for (FTPFile bFile : subFiles) {
          // System.out.println(subDir + "\t" + bFile.getName());
          map.put(bFile.getName(), collectionName + "/");
        }
      } else {
        // System.out.println(rootDir + "\t" + aFile.getName());
        map.put(aFile.getName(), "");
      }
    }
    return map;
  } // makeMapName2Subdir

  protected boolean determineFtpSource(String host, FTPClient ftp, EnsemblGenome genome, int retryCount) throws IOException, InterruptedException {
    for (int attempt = 1; attempt <= retryCount; attempt++) {
      try {
        // this could happen if we're retrying
        if (!ftp.isConnected()) {
          // try to recreate the connection
          if (!loginFtp(ftp, host)) throw new IOException("Unable to log in to " + host);
        }
        return determineFtpSource(ftp, genome);
      } catch (IOException e) {
        logger.log(Level.WARNING, "Failed to determine file for " + genome + " " + attempt + " of " + retryCount, e);
        if (attempt == retryCount) throw e;
        if (ftp.isConnected()) {
          try {
            ftp.disconnect();
          } catch (IOException e2) { /* ignore */ }
        }
        Thread.sleep(TimeUnit.SECONDS.toMillis(10));
      }
    }
    return false;
  }

  protected boolean determineFtpSource(FTPClient ftp, EnsemblGenome genome) throws IOException {
    progress.setTask("Determining full URL for " + genome, 0, -1);
    PatternFileFilter filter = new PatternFileFilter(genome.getFilePattern());
    logger.log(Level.INFO, "Looking for " + genome + " sequences at ftp://" + genome.getRemoteHost() + genome.getRemoteDir());
    FTPFile sequenceFile = filter.best(ftp.listFiles(genome.getRemoteDir(), filter));
    if (sequenceFile == null) {
      logger.log(Level.WARNING, "Skipping ftp://" + genome.getRemoteHost() + " " + genome + " as no sequence found.");
      return false;
    }
    logger.log(Level.INFO, "Using file " + sequenceFile.getName() + " for " + genome);
    genome.setRemoteInfo(sequenceFile.getName(), sequenceFile.getSize());
    return true;
  } // determineFtpSource

  /**
   * Checks if a entry in the database exists for the combination of
   * category name, listing name, sequence alphabet and sequence edition. Can also
   * check for editions that are newer or within delta of the given edition.
   * In some rare cases a sequence of interest may have been marked obsolete, for
   * example if the latest edition has been manually deleted, this also allows
   * removing the obsolete flag.
   * @param genome contains the information about the genome to be queried.
   * @param deobsolete should an existing genome matching the criteria be deobsoleted?
   * @param exact does the match have to be exact or are newer editions allowed?
   * @param delta the allowance made for editions to be older.
   * @return true if the database contains a sequence matching the criteria.
   * @throws SQLException if the database access doesn't work.
   */
  protected boolean sourceExists(Source genome, boolean deobsolete, boolean exact, long delta) throws SQLException {
    boolean exists = false;
    Long fileId = null;
    Boolean obsolete = null;
    Connection conn = null;
    dbLock.readLock().lock();
    try {
      conn = dataSource.getConnection();
      PreparedStatement ps = conn.prepareStatement(
          exact ? SELECT_SEQUENCE_FILE_BY_INFO : SELECT_SEQUENCE_FILE_BY_INFO_NEWER);
      ps.setString(1, genome.getCategoryName());
      ps.setString(2, genome.getListingName());
      ps.setInt(3, genome.guessAlphabet().toInt());
      ps.setLong(4, genome.getSequenceEdition() - delta);
      ResultSet resultSet = ps.executeQuery();
      if (resultSet.next()) {
        fileId = resultSet.getLong(3);
        obsolete = resultSet.getInt(4) != 0;
        exists = true;
      }
      resultSet.close();
      ps.close();
      if (deobsolete && exists && obsolete) {
        dbLock.writeLock().lock();
        try {
          ps = conn.prepareStatement(MARK_SEQUENCE_FILE);
          ps.setInt(1, 0); // set obsolete = 0
          ps.setLong(2, fileId);
          ps.executeUpdate();
          ps.close();
        } finally {
           dbLock.writeLock().unlock();
        }
      }
      conn.close();
      conn = null;
    } finally {
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) {
          logger.log(Level.SEVERE, "Unable to close database connection after error.");
        }
      }
      dbLock.readLock().unlock();
    }
    return exists;
  }

  protected boolean sourceExists(Source genome, boolean deobsolete) throws SQLException {
    return sourceExists(genome, deobsolete, true, 0);
  }

  protected boolean loginFtp(FTPClient ftp, String host) throws IOException {
    int reply;
    ftp.setConnectTimeout(CONNECT_TIMEOUT);
    ftp.setDefaultTimeout(REPLY_TIMEOUT);
    ftp.setAutodetectUTF8(true);
    logger.log(Level.INFO, "Attempting connection (" + host + ").");
    progress.setTask("Attempting connection to " + host, 0, -1);
    ftp.connect(host);
    reply = ftp.getReplyCode();
    if (!FTPReply.isPositiveCompletion(reply)) {
      logger.log(Level.SEVERE, "Connection attempt failed (" + host + ").");
      ftp.disconnect();
      return false;
    }
    logger.log(Level.INFO, "Logging in (" + host + ").");
    progress.setTask("Logging in", 0, -1);
    if (!ftp.login("anonymous", "meme-suite@uw.edu")) {
      logger.log(Level.SEVERE, "Log in attempt failed (" + host + ").");
      ftp.disconnect();
      return false;
    }
    logger.log(Level.INFO, "Setting file type to binary (" + host + ").");
    progress.setTask("Setting file type to binary", 0, -1);
    if (!ftp.setFileType(FTP.BINARY_FILE_TYPE)) {
      logger.log(Level.SEVERE, "Unable to set file type to binary (" + host + ").");
      ftp.disconnect();
      return false;
    }
    ftp.enterLocalPassiveMode();
    return true;
  }

  protected boolean downloadFtpSource(FTPClient ftp, FtpSource source, boolean supportResume, int retryCount) throws IOException, ExecutionException, InterruptedException {
    for (int attempt = 1; attempt <= retryCount; attempt++) {
      try {
        // this could happen if we're retrying
        if (!ftp.isConnected()) {
          // try to recreate the connection
          if (!loginFtp(ftp, source.getRemoteHost())) throw new IOException("Unable to login to " + source.getRemoteHost());
        }
        return downloadFtpSource(ftp, source, supportResume);
      } catch (IOException e) {
        logger.log(Level.WARNING, "Failed to download source " + source + " attempt " + attempt + " of " + retryCount, e);
        if (attempt == retryCount) throw e;
        if (ftp.isConnected()) {
          try {
            ftp.disconnect();
          } catch (IOException e2) { /* ignore */ }
        }
        Thread.sleep(TimeUnit.SECONDS.toMillis(10));
      }
    }
    return false;
  }

  protected boolean downloadFtpSource(FTPClient ftp, FtpSource source, boolean supportResume) throws IOException, ExecutionException, InterruptedException {
    checkWorkerTasks();
    progress.setTask("Downloading " + source, 0, source.getRemoteSize());
    // create the downloads directory if it does not exist
    File downloadDir = new File(dbDir, "downloads");
    if (downloadDir.exists()) {
      if (!downloadDir.isDirectory()) {
        throw new IOException("Unable to create download directory \"" +
            downloadDir + "\" as a file with that name already exists!");
      }
    } else if (!downloadDir.mkdirs()) {
      throw new IOException("Unable to create download directory \"" +
          downloadDir + "\" as the mkdirs command failed!");
    }
    File target = new File(downloadDir, source.getLocalName());
    logger.log(Level.INFO, "Downloading " + source + " from " +
        source.getRemoteUrl() + " to \"" + target + "\".");

    long lastCommandTime;
    FileOutputStream out = null;
    InputStream in = null;
    boolean reconnectAfterTransfer = false;
    try {
      long downloaded = (target.exists() && supportResume ? target.length() : 0);
      progress.setTaskPartDone(downloaded);
      out = new FileOutputStream(target, supportResume);
      if (supportResume) ftp.setRestartOffset(downloaded);
      in = ftp.retrieveFileStream(source.getRemotePath());
      lastCommandTime = System.currentTimeMillis();
      FtpCopier copier = new FtpCopier(in, out);
      copier.start();
      while (copier.isAlive()) {
        if (System.currentTimeMillis() - lastCommandTime > KEEP_ALIVE_TIMEOUT) {
          lastCommandTime = System.currentTimeMillis();
          try {
            logger.log(Level.FINE, "Sending Keep-alive NOOP to " + source.getRemoteHost());
            ftp.sendNoOp();
            logger.log(Level.FINE, "Keep-alive NOOP to " + source.getRemoteHost() +
                " acknowledged.");
          } catch (SocketTimeoutException e) {
            reconnectAfterTransfer = true;
            logger.log(Level.FINE, "Keep-alive NOOP to " + source.getRemoteHost() +
                " ignored!", e);
          }
        }
        checkWorkerTasks();
        copier.join(TimeUnit.SECONDS.toMillis(1));
      }
      copier.check();
      out.close();
      out = null;
      in.close();
      in = null;
    } finally {
      if (in != null) {
        try {
          in.close();
        } catch (IOException e) {
          logger.log(Level.WARNING, "Ignoring failure to close FTP InputStream.", e);
        }
      }
      if (out != null) {
        try {
          out.close();
        } catch (IOException e) {
          logger.log(Level.WARNING, "Ignoring failure to close downloaded file OutputStream.", e);
        }
      }
    }
    boolean isReplying = false;
    boolean transferOk = true; // assume true until told otherwise
    try {
      transferOk = ftp.completePendingCommand();
      isReplying = true;
    } catch (SocketTimeoutException e) { /* ignore */
      reconnectAfterTransfer = true;
      logger.log(Level.FINE, "Handling Socket Timeout from completePendingCommand", e);
    }
    if (reconnectAfterTransfer) {
      // Unfortunately our command channel will probably have been closed...
      logger.log(Level.INFO, "Disconnecting to recreate closed control channel");
      if (isReplying) ftp.logout();
      if (ftp.isConnected()) {
        try {
          ftp.disconnect();
        } catch (IOException e) {
          logger.log(Level.FINE, "Disconnecting caused exception", e);
        }
      }
      // try to recreate the connection
      loginFtp(ftp, source.getRemoteHost());
    }
    // check for success, if we somehow kept the connection
    if (!transferOk) {
      logger.log(Level.WARNING, "Transfer command failed (" + source + ").");
      return false;
    }
    // check the file size
    if (target.length() != source.getRemoteSize()) {
      logger.log(Level.WARNING, "Transfered file size does not match remote (" + source + ").");
      return false;
    }
    logger.log(Level.INFO, "Download of " + source + " complete.");
    source.setSourceFile(target);
    return true;
  }

  protected boolean downloadHttpSource(URL url, FtpSource source) throws IOException, ExecutionException, InterruptedException {
    progress.setTask("Downloading " + source, 0, source.getRemoteSize());
    // create the downloads directory if it does not exist
    File downloadDir = new File(dbDir, "downloads");
    if (downloadDir.exists()) {
      if (!downloadDir.isDirectory()) {
        throw new IOException("Unable to create download directory \"" +
            downloadDir + "\" as a file with that name already exists!");
      }
    } else if (!downloadDir.mkdirs()) {
      throw new IOException("Unable to create download directory \"" +
          downloadDir + "\" as the mkdirs command failed!");
    }
    File target = new File(downloadDir, source.getLocalName());
    logger.log(Level.INFO, "Downloading " + source + " from " +
        source.getRemoteUrl() + " to \"" + target + "\".");

    InputStream in = null;
    FileOutputStream out = null;
    try {
      // Copy the file to the downloads directory.
      URLConnection connection = url.openConnection();
      HttpURLConnection httpConnection = (HttpURLConnection)connection;
      httpConnection.setRequestProperty("Content-Type", "application/octet-stream");
      InputStream response = connection.getInputStream();
      int responseCode = httpConnection.getResponseCode();
      if (responseCode != 200) {
	throw new RuntimeException("Response code was not 200. Detected response was " + responseCode);
      }
      //char[] buffer = new char[8192];
      byte[] buffer = new byte[2048];
      int read;
      //in = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      //in = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      in = response;
      out = new FileOutputStream(target, false);
      while ((read = in.read(buffer, 0, buffer.length)) > 0) {
	out.write(buffer, 0, read);
	progress.addTaskPartDone(read);
      }
    } finally {
      if (in != null) {
        try {
          in.close();
        } catch (IOException e) {
          logger.log(Level.WARNING, "Ignoring failure to close HTTP InputStream.", e);
        }
      }
      if (out != null) {
        try {
          out.close();
        } catch (IOException e) {
          logger.log(Level.WARNING, "Ignoring failure to close downloaded file OutputStream.", e);
        }
      }
    }
    boolean transferOk = true; // assume true until told otherwise
    // check for success, if we somehow kept the connection
    if (!transferOk) {
      logger.log(Level.WARNING, "Transfer command failed (" + source + ").");
      return false;
    }
    // check the file size
    if (target.length() != source.getRemoteSize()) {
      logger.log(Level.WARNING, "Transfered file size does not match remote (" + source + ").");
      return false;
    }
    logger.log(Level.INFO, "Download of " + source + " complete.");
    source.setSourceFile(target);
    return true;
  }

  public synchronized void checkWorkerTasks() throws ExecutionException, InterruptedException {
    if (Thread.currentThread().isInterrupted()) throw new InterruptedException();
    Future task;
    while ((task = workerEcs.poll()) != null) {
      // tasks should be returned in roughly queue order so this operation shouldn't be as time intensive as it looks...
      if (!workerTasks.remove(task)) {
        throw new Error("An task was returned that we don't own!");
      }
      try {
        task.get();
      } catch (CancellationException e) {
        // ignore
      }
    }
  }

  public synchronized void enqueueSequences(SequenceProcessor sequenceProcessor) {
    Future<Void> future =  workerEcs.submit(sequenceProcessor);
    workerTasks.add(future);
  }

  public synchronized void cancelWorkerTasks() {
    ListIterator<Future<Void>> iterator = workerTasks.listIterator(workerTasks.size());
    while (iterator.hasPrevious()) iterator.previous().cancel(true);
  }

  public synchronized void waitForWorkerTasks() throws ExecutionException, InterruptedException {
    if (Thread.currentThread().isInterrupted()) throw new InterruptedException();
    Future task;
    while (!workerTasks.isEmpty() && (task = workerEcs.take()) != null) {
      // tasks should be returned in roughly queue order so this operation shouldn't be as time intensive as it looks...
      if (!workerTasks.remove(task)) {
        throw new Error("A task was returned that we don't own!");
      }
      try {
        task.get();
      } catch (CancellationException e) {
        // ignore
      }
    }
  }

  protected static Properties loadConf(Class<?> owner, File dbDir, String name) {
    Properties conf = new Properties();
    File confFile = new File(dbDir, "conf" + File.separator + name);
    try {
      FileUtils.copyInputStreamToFile(owner.getResourceAsStream(name), confFile);
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Failed to copy configuration file from defaults", e);
      return conf;
    }
    try {
      conf.load(new FileInputStream(confFile));
    } catch (IOException e) {
      logger.log(Level.SEVERE, "Failed to load configuration file", e);
    }
    return conf;
  }

  protected static Integer getConfInt(Properties conf, String key, Integer min, Integer max, Integer defVal) {
    String strVal = conf.getProperty(key);
    not_valid:
    if (strVal != null) {
      Pattern INT_RE = Pattern.compile("^\\s*([+-]?\\d+)\\s*$");
      Matcher m = INT_RE.matcher(strVal);
      if (!m.matches()) {
        logger.log(Level.WARNING, "Value for " + key + " was not an integer:\n" + strVal);
        break not_valid;
      }
      int val;
      try {
        val = Integer.parseInt(m.group(1), 10);
      } catch (NumberFormatException e) {
        logger.log(Level.WARNING, "Value for " + key + " was not parsable as an integer", e);
        break not_valid;
      }
      if (min != null && val < min) {
        logger.log(Level.WARNING, "Value for " + key + " was too small " + val + " < " + min);
        break not_valid;
      }
      if (max != null && val > max) {
        logger.log(Level.WARNING, "Value for " + key + " was too large " + val + " > " + max);
        break not_valid;
      }
      return val;
    }
    return defVal;
  }

  private class FtpCopier extends Thread {
    private InputStream in;
    private OutputStream out;
    private Throwable throwable;
    public FtpCopier(InputStream in, OutputStream out) {
      this.in = in;
      this.out = out;
    }
    @Override
    public void run() {
      int read;
      byte[] buffer = new byte[2048];
      try {
        while ((read = in.read(buffer)) != -1) {
          out.write(buffer, 0, read);
          progress.addTaskPartDone(read);
        }
      } catch (Throwable e) {
        this.throwable = e;
      }
    }

    public void check() throws IOException {
      if (throwable != null) {
        if (throwable instanceof IOException) {
          throw (IOException)throwable;
        } else if (throwable instanceof RuntimeException) {
          throw (RuntimeException)throwable;
        } else if (throwable instanceof Error) {
          throw (Error)throwable;
        } else {
          throw new Error(throwable);
        }
      }
    }
  }

}
