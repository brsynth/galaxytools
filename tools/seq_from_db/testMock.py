import psycopg2
import subprocess
import time


def docker_container_exists(name):
    """Vérifie si le conteneur Docker existe (arrêté ou démarré)."""
    result = subprocess.run(
        f"docker ps -a -q -f name={name}",
        shell=True, capture_output=True, text=True
    )
    return bool(result.stdout.strip())


def docker_container_is_running(name):
    """Vérifie si le conteneur Docker est en cours d'exécution."""
    result = subprocess.run(
        f"docker ps -q -f name={name}",
        shell=True, capture_output=True, text=True
    )
    return bool(result.stdout.strip())


def docker_start_container(name):
    """Démarre un conteneur Docker existant."""
    subprocess.run(f"docker start {name}", shell=True, check=True)


def docker_run_container(name, password="RK17", port=5432):
    """Crée et démarre un nouveau conteneur PostgreSQL."""
    subprocess.run(
        f"docker run --name {name} -e POSTGRES_PASSWORD={password} -p {port}:5432 -d postgres",
        shell=True, check=True
    )


def docker_stop_container(name):
    """Arrête un conteneur Docker."""
    subprocess.run(f"docker stop {name}", shell=True, check=True)


def wait_postgres_ready(host="localhost", port=5432, timeout=30):
    """Attend que Postgres accepte les connexions."""
    import socket
    start = time.time()
    while time.time() - start < timeout:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            try:
                sock.connect((host, port))
                return True
            except ConnectionRefusedError:
                time.sleep(1)
    raise RuntimeError("PostgreSQL is not ready after timeout")


def create_db_and_insert_data():
    container_name = "test_fragments_db"

    # 1. Vérifier si conteneur existe
    if docker_container_exists(container_name):
        print(f"Conteneur Docker '{container_name}' existe.")
        # 2. Démarrer s'il n'est pas déjà lancé
        if not docker_container_is_running(container_name):
            print(f"Démarrage du conteneur '{container_name}'...")
            docker_start_container(container_name)
        else:
            print(f"Conteneur '{container_name}' déjà en cours d'exécution.")
    else:
        print(f"Création et démarrage du conteneur '{container_name}'...")
        docker_run_container(container_name)

    # 3. Attendre que PostgreSQL soit prêt
    wait_postgres_ready()

    # 4. Connexion et création DB + table + insertion
    conn = psycopg2.connect(
        dbname='postgres', user='postgres', password='RK17',
        host='localhost', port='5432'
    )
    conn.autocommit = True
    cursor = conn.cursor()

    cursor.execute("SELECT 1 FROM pg_catalog.pg_database WHERE datname = 'test_fragments_db';")
    exists = cursor.fetchone()

    if exists:
        print("Base de données 'test_fragments_db' déjà existante, suppression...")
        cursor.execute("DROP DATABASE test_fragments_db;")

    cursor.execute("CREATE DATABASE test_fragments_db;")
    print("Base de données 'test_fragments_db' créée.")

    cursor.close()
    conn.close()
    # Connect to the PostgreSQL container (default 'postgres' database for setup)
    conn = psycopg2.connect(
        dbname='postgres',  # Default database
        user='postgres',  # Default user
        password='RK17',  # Password from Docker environment
        host='localhost',  # Running locally on the default Docker network
        port='5432'  # Default PostgreSQL port
    )

    conn.autocommit = True  # Necessary to create a database
    cursor = conn.cursor()

    # Check if the test database already exists
    cursor.execute("SELECT 1 FROM pg_catalog.pg_database WHERE datname = 'test_fragments_db';")
    exists = cursor.fetchone()
    
    if exists:
        print("Database 'test_fragments_db' already exists, dropping it...")
        cursor.execute("DROP DATABASE test_fragments_db;")

    # Create the new database for testing
    cursor.execute('CREATE DATABASE test_fragments_db;')
    print("Database 'test_fragments_db' created.")

    cursor.close()
    conn.close()

    # Now connect to the new test database
    conn = psycopg2.connect(
        dbname='test_fragments_db',
        user='postgres',
        password='RK17',
        host='localhost',
        port='5432'
    )

    cursor = conn.cursor()

    # Create the 'sample' table instead of 'fragments'
    cursor.execute('''
        CREATE TABLE sample (
            fragment TEXT PRIMARY KEY,
            sequence TEXT,
            annotation TEXT,
            metadata_1 TEXT,
            metadata_2 TEXT
        );
    ''')

    # Insert mock data
    cursor.executemany('''
    INSERT INTO sample (fragment, sequence, annotation, metadata_1, metadata_2)
    VALUES (%s, %s, %s, %s, %s);
''', [
    (
        'part_A',
        '''ORIGIN
        1 tcacaggatg gtccaacgaa actaggcttt agacgaggga tgaatgaccg acccccactc
       61 gtggcactaa cggacagact tccctgacgg ttattcgacc attaaagtca gacatgcggg
      121 ggtgaataaa ttagccaaat tgtgtcgaag aaaagacgtg cggctggcac ataaggcagt
      181 cttgatccta gtcttgcagg gatgcacgta agtcgcctca attaactgca gccgagctcc
      241 aggttaccaa agaccctagt atgccagggc ctaacggttg gagtatatta tgggtacgca
      301 atagtgcgga agttaacctg ggcaacatcc aggtgagagg ttggacggaa gcgacagtaa
      361 gtggccatag actgccgagt cgtgttaatg aatcgctata cgcccatgga gttgtggggt
      421 cgttttatcc gagtaggggc ccgctgacta cttcgtccag acaatatgcc gtcttcaata
      481 gtctacctga gagtcatgcc ggcatttccg acgctgagtg aaacccgcgt agccaggcga
      541 aatttgcatc ttgaaatacc actgcagatc agccagtaag gcccatataa gggcgctggt
      601 gttctggcga cagataaagt gttatctaat gtaacccgcg gacttttaga ggatacttga
      661 atgcgggcat atcgtcccac cagcgtcacg tggtcgggct agcggcagac aaactctgcc
      721 gatgtttttc tttgccgcga tagcgggcta gtcaattcac tatccggcga tgaagcatag
      781 attgttatcg cgcttatgcg gaggataatc aagtcttggc agaacctgct cgctcatgca
      841 ctggtgaggc gtggttatcc taacaaccgc ctaggacgaa gaatgggctc ggtagggaca
      901 gttcgtgcgt ttagcttcgt cccatcttaa tgctgtctgg agggaggcta catgcgaaca
      961 gaagccgtgg gcgaaagttg ttgatccggg tctaacacgg ataagggcca tggtgaacgc
     1021 atcatgcgat acttcttggt agttcttcta acatgaccgt gccatagccc ttaccctttg
     1081 aagaagttaa cctatctgcc gtctcctgca cgaatagcga ctggactcgg gtattttgga
     1141 tacctcacga aagcacttcg atcggccgaa gcatggatat ttcgccggac gggcccgaat
     1201 agagagcctt gtatgtcggt tagtacaacg cagtcctgga gacatctacg cggatgggcc
     1261 taggggctgg acttaacatt gggtaacgta cctggtccaa agtgaatgca aagcactttt
     1321 acaaagcgcg ggtccccgtg ggtgtttagg gtagaagatt ggcggatgct acgaacgatc
     1381 ccgctttgaa actatcatta catcgtgtaa aagagacact taacaaaggc caataaactg
     1441 ccagtaagaa tcgctcagtg cggtgctggg gacgctaagt aggggcaaca gccagtgaag
     1501 gacgtgaccg acctttctca gataagatat gctggcgtct atctaataag catagtgaaa
     1561 aaccaaccat ttcacttaca cgaagtacat ttgcattgct agtaaagacg cctaaacaga
     1621 agtgcccttg catgctgtat gtctatagtc cttagggaag catcagccct tctacttatt
     1681 cgaggtctga gaaaccctgg acaagctccg aattattcaa tgtgcctgtc tccgaggtta
     1741 gatagcgcta tgctcttaag agttgcacag aatgaccatc ctggaatgtc cctggagggg
     1801 tctaggtatg ctgatcgaag ggtgctctaa ggacttgacg tgcgtccgag gagggtgctg
     1861 cctccttcgc ctttagatcc aacgcggatc acatgcgcgt gagctaatag gatcaccttc
     1921 tgctccgatt tttaccctcc tgggtcactt ccgaatgagg tagcgggcga aatataatgt
     1981 ctccactcgt aggtgttccc tgtatgtgaa gctctatagt ggacaaaggt ttgatgaact
     2041 agcccccgta tacgctctca ccgacggacg cgggggtctg ttatttgaag catcatacat
     2101 gcgaaggtgc cttctcagca acgaaaggta gtgggagtgt acaagttcaa tgcgccgcca
     2161 taggtctgag tatacaaggg gatgccccca tccacaacgg gattggctac ccggagagct
     2221 ggctccgctc caacaaataa ttatattaac ctattggaat tccacctgca tatcagagga
     2281 gagagacctt tacggctatt ctgtttaccg gatccatcgg taccaaggat cagaaagtga
     2341 cacggtttga acgggttgtt gtaatacttt gagtatacct ctgacgctga gcgtgtcgtc
     2401 ctgagcgcag actcaataac atagcagtcc gacatcgccg tgatatgtaa atgcaacgaa
     2461 tttaggtctt gactcggtct accatgtcaa aagggtagcc agatttcagc gcgaaattga
     2521 actttgtgtt tagtgtgggg tcctcggtta caaaatagga tcagacatgt gtgattttgg
     2581 taacctagtc tggcagtccg acagacttcg ctatgatttg atggggccgg cctataattg
     2641 gcttgcgcaa cccgctcatc tcgggcgtgt tttacttcct gcggtcccca cgccctattt
     2701 tcgggccagc tgtaggtgct agagtgaatg ctggcgaata agattccccg ctatttggcg
     2761 cctcgccaca gctctggcac tatgggggga gtttctctgt tccttaaaca gcacccgttt
     2821 ttgaggtgta ttggtttcgg ttctgcatta ggcaattcgt accgtacaat caattacgac
     2881 acattggcgg cagttatcag ctacccatcg caaagcacac acccacatgt atctattttt
     2941 cgcaaattcc aaaagcttcg attgagattg catcggtagt ccctcagaca tgtcgtaatc
     3001 gaatgcctct tgttccatga gagagagaag tatggcgcga accgctctgc ctttaatttg
     3061 gttctaccat ccacgagttt aaggggcata accctgccca gcactttccg aggctcacgt
     3121 tcatgctacg gtagcacctt tctgcgggtc tcacgctgtc aatatgcagg tgctgcagga
     3181 atttgtctcc aatcgacttg agatatcgca agcatgaaat tatattagac acgccagaga
     3241 acttgggaag cagcactggt agtgatagca acccgagtac agtaacgagt gagcttctga
     3301 tcatgagctc tcctacggcg tcaatgcgac gaatgcccag catgcactct cgctatccat
     3361 gcctgctagg gtggcattat gctcaggaac agttgtagct tggatatcgt ctagatgaaa
     3421 tacctggaca ctggttagcg tcgtcaagca ccaaggacat tcacacgctc gcggtctttc
     3481 gtctccctaa gcgttcggca gtcgggcgtg aagaggttgt aatcagacgg aacaaagcct
     3541 gaaaaaattc cagcgacgta gtattcatga tcctgtacca tctgtagccg ccgcggcgca
     3601 cgattgaatg taggctacta accccatccg tgttagcgat gtgagtttct accgcaacga
     3661 atgctcaagc gaaccttctt ctttcgtccg caacccacaa gccgtggtta tgacagctaa
     3721 attgtcccag acatcccttt attacacaag agctccagcg gaatacctag tcacagcggt
     3781 aatgacacaa agctcttagt tagtccaggg actacttctg tctacagcac atcacactca
     3841 ttatcagcat cagtgtagag acggagaaca tgggctatcc tataccaaga tccgccatct
     3901 aaacatttga agtttcccgt cttctataac ttagcactcg acgctattct gctgagtgcg
     3961 cttagtctgt agcgacttgc gaaatccata aactgagaat tgaaagagag tgcataaccg
     4021 aaacctttgt ggcatatttc cgttgaaacg taccagaaca gccgtttagt gcggaacata
     4081 cagtatcctg ataaagcact caacccaaca gacaccctat gccgatagcg ggatgctaac
     4141 aagtatagtc atgatgattt ctcggacagc ggggtttggt acagctgcaa tccgtgattt
     4201 aaattcggac ctctgcacac accgacggtt acccatatcc tctacggctg taggaagttt
     4261 taccttggat gtcattctga tttcggcgta tc
''',
        '''LOCUS       part_A                  4292 bp    ds-DNA  circular UNK 28-OCT-2019
DEFINITION  .
ACCESSION   part_A
VERSION     part_A
KEYWORDS    "creator:SynthSys Center".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     RBS             2332..2343
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2208..2238
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     promoter        2289..2323
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#85dae9"
     misc_feature    3162..3178
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -35_signal      228..233
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     gene            1009..2028
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -35_signal      51..56
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             2350..3027
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff00ff"
     CDS             complement(3376..4191)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#993366"
     misc_RNA        87..639
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2257
                     /label="another feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     terminator      3119..3146
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    2053..2100
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    2278..2281
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_RNA        complement(90..197)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    726..734
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2274..2277
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
     rep_origin      51..639
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2161..2185
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -10_signal      206..211
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    complement(3292..3332)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     rep_origin      760..1008
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     terminator      3032..3103
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    3158..3161
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    3173..3178
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     gene            complement(3376..4191)
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2258..2262
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     misc_feature    3342..3365
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    3204..3227
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2256..2273
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -10_signal      72..77
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2029..2052
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             1009..2028
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#993366"
     misc_feature    3237..3265
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    3154..3157
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_feature    3180..3290
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    2212..2236
                     /label="feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    1..24
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
''',
        'Metadata1 for part_A', 'Metadata2 for part_A'
    ),
    (
        'part_B',
        '''ORIGIN
        1 cgggtctagc cacgcggctg aaaatgcgcc tcaccaagtt cagagagttg taatttaaag
       61 gaccaaagaa cgtccattat gaagacacca atgagaaaag ttcattgcga tcgcagctta
      121 gactgacgtg aagttccaac atttggcata cgtacttgtt aagttcctcc ggctatatgc
      181 gctaggcgca ggtatattag ccggacatct gcaaggttga cgctgaagcc ttatcggtaa
      241 tgcaagcaca acggtattgc agacacacaa taggaagtta tgttaacagg agagggtgat
      301 gccgatgatt acgccccgcg tcgacagtac aataaaacag ctagatcatg tcgaaagtgt
      361 aggcaaaaca ggtaggtgca ggacaggagg gttccctcaa gtcggcttgc tttatgtggc
      421 gaatctgata cccataccgg ctatacatcg gccagggaat atcacagatt ttaaatgagc
      481 ggagatagac ctatgggacc gcatcgatgc aagacatcgc ctacagatac ataatctgga
      541 agacgcagat atcagtacgc tccttcctgg gaggcgtggt ccttggcacc agtcggccac
      601 cccttagcgt ttcacctaag caacgaatcc aggactcctc ttttagcgtt tgggtgtcag
      661 gcgcatcagg agtggttaaa tgaacgctgt acttcgcaag ataattgtga agatctctgc
      721 tattccttaa cattaccacc aactttgcct acatatggcg taatcaaaaa tgcagcgggg
      781 tgctcgatcg tgaaccctca atcacgatgt gctatggcgt agccacacaa catagttccg
      841 gtctggaatt cggtctctgg agagcatacg gtgctctggg ctcacacgta cttcctcttg
      901 tcacgcgcct ttgtagagca ttagggcgga ggtcactggg acggggatca gagagcttaa
      961 aatcccaagt agcagtgcct aaggtgtaga atcttgtgca gatcttcgtc tggcccgaaa
     1021 cctgggcctc cacaacaatg tcgggtcttt tataatttga aaacatcgcg ccttacccgc
     1081 accgcgcagg acgattagac acagggggtc atgtccatgg ttcggggtaa gccaaattta
     1141 catcgagtgt gggcttcatt gaggagaggg cgttacgata cccgttaaat gaaaacccgg
     1201 acagtctcta actgataggc ctgataatcg ccaaatgaag tcgcaattag cagtacatga
     1261 ctggcagata tgtcgccgga gacggtatgg gggacctatc ctaggttggc acgcgagcgc
     1321 gaatgacgcc tatcttatcc ccctgcacgt gatggtatcg tcccaactat cttttacctc
     1381 atcaacgttt tcctgcataa tctgcaagat ggagagtccg cgagtgcggg gagtacgaga
     1441 tcacaaatct gccctgaaga caaagctcgc aggccgtact caggaaatag catcgactca
     1501 cacatatcgg gaacggcaaa gagtgcatcg gccaatgcaa actattgaac gtccatgaaa
     1561 ttgtgtcaaa ctactagtgt gtggacttgc gtcaaggcca tggcgcaata ctcaaagccc
     1621 cggtatttat cccccgcaga gtatccgtgc ggatggagtt ccacagctct ttggacatgc
     1681 aatgcgtggc gctacatgct ccgctgatac agaatatcta atgtctcgcc tcttgccgtc
     1741 gaccctacaa ctactgcacg ggcagagtgg acgcgacccc tgagcggtgc agaaccttag
     1801 tataaggcaa gggtggcacg gccgagccgc acttgctcct actaggccgt gctcgaagct
     1861 attgatgcaa gcacagtatt cccgatttgt aaagcaaggt gcccgtagcc atacacccag
     1921 ctactcaccg gccaatttcg accctctcta tctcagacct ccgcatattc ccaaactcct
     1981 gctcaagtat gaaatgaccg gctcctccta aggcctgcgg acatggagtc atgcgtggaa
     2041 tgccgtgtag gtagggtaca cacagcctta ccatgtaggg caataaatta tagactactg
     2101 acattgattc cttgagggcc cgttcacttg tgggggaacc aggaaggcgc cgatagtcat
     2161 cacttaaaca catgcgacgc atccatttca agaaagacac gctacgggtc tacaccgcct
     2221 atagttctaa tttgggcccc ggtacgcaca tcgatgtact ctccattaca gcatagcaca
     2281 tcagcggccc tagttagcgc ttctttttgt tccaacattc aatagcgggg actgcattat
     2341 cctcacccta tttgttgcat aggtatttgg ccgtatactt gggaaagagc aaccacggga
     2401 gtgcactaga aactgattca cccatgttat aacgcgttgc tcgcaagtcg tactgtggcc
     2461 acgcacaggt caggcaatcg caaccatgtc ttagcagttt gtaagaagac ctttaaacat
     2521 gtacctcaca gcaatctcgg aacttttaat tttcctgagt caatgacatg ctagctccgt
     2581 ccacttgtaa tcaccgtcat aacacagtta gcctgcaaaa gggcgtgtac accttgccac
     2641 tttaaagctg accgaataac ctacctacag gagatgctta gtgccacgga atgatctcgc
     2701 aatgcttagt cttacaccca gtgaatcaac ccaagacctc gccgactgaa cggtattgag
     2761 tccaatcccg agaggcaatt tatggccgtt agtgacagtc ggtcgctgat agaataaggc
     2821 agactcgtct taagcgggct gctttaagag aaaatggagc cctgacgaag aagtggatga
     2881 accaggtttt gacaacgtcc catattggta atgtaggctc gttatcaaga gataggttta
     2941 gaaggctcaa tgctcgacgt aaaagaagtc ttaagctcta tgagccagtg atgagagcta
     3001 attatagggc gcattcgttg gggtaaggtt gcgaattgtc atttaactaa aaccgacggc
     3061 aatgtgagac cctgcagtcc gcctacacct ggctcactga ggcggcaacc gttatcgaag
     3121 cagcatgtat tggtagtacc tgtcgccgat ttccttggga gtgcgctgtt gcgataacgc
     3181 ttatttatgg tttgcgcggg atgtatctaa tgcggccgcc gcccgatcag cctgaacact
     3241 gcggtcatga gccggtccat catgcagttc gcctgtaatt ataacctatt tgttaatggg
     3301 gcttacagag gagcttccaa cgcatctacg agtaacggct gggatcacaa cttagccacc
     3361 ttgtagcagg gacataccgc tttgcgacta ctgcgattct ctaacgagtc ttctgattgg
     3421 tcctatgctc tttcaaaatc tcgataccgt agccaaccct gctgcccacc ttaaactatg
     3481 acaggacaca ccattgtgcc aacctacacc ttaccagtag ttctgctctc caccactaat
     3541 acgatgtcag aggccagaag gagcaccgtc aacaaacgcg ataatcaggg cggcgaagat
     3601 ggaagtagct tgggttgtcg gtgggcaggt caggcgttag tgttatttga ataccaactg
     3661 ttagaccaca gcgactcgaa cagtccgcaa gtcaagtgcg gtcggctttt acaggtggcg
     3721 gattttcatt ccaccattcg ggtctaggtg tggtactaac cgcggtcata cgtggagttg
     3781 aaccgcctga cgcactaagc tatgcgggca acaatatcgc caaagtcatt gtgttgggtt
     3841 ttcacaacag tcacttcgct cacggccgaa ttaaatatct ctaactatta agccattttc
     3901 caaccatcgg ctatcccagc ggttcaatcc gaatcgttac cagtgggtag ctgcagagga
     3961 tcctaagggc gggtcccgct tactgtaact cgtcatgtct cgtggtctac caattcggtt
     4021 agtattcgtg ccatgtatct cataagggtt gcagtttcgt cattgtcgta gttccttaat
     4081 ggactggtac gacccggttg ttatacaacc atcaatcgtt aacggacatt agtgagccga
     4141 ggagaagacc catcagcgga taacccgcca ttatcaccga tactgacata tgagtgccat
     4201 actgggactg gcttcacgtg agggacctca aagtgagtcg tacaaggcag ctcctcaaa
''',
        '''LOCUS       part_B                  4259 bp    ds-DNA  circular UNK 29-OCT-2019
DEFINITION  .
ACCESSION   part_B
VERSION     part_B
KEYWORDS    "Source" "Sequence" "Freezer" "Shelf" "Box" "Concentration" "Use".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     primer_bind     531..550
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     terminator      3078..3135
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     CDS             complement(join(4216..4259,1..616))
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#f58a5e"
     exon            2863..3060
                     /label="feature"
                     /ApEinfo_fwdcolor="#75c6a9"
     gene            2863..3060
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
     terminator      4101..4195
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#9eafd2"
     primer_bind     complement(3410..3429)
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     primer_bind     778..796
                     /label="another feature"
                     /ApEinfo_fwdcolor="#75c6a9"
     misc_feature    859..862
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     primer_bind     complement(3159..3176)
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    863..3060
                     /label="another feature"
                     /ApEinfo_fwdcolor="#f58a5e"
     rep_origin      complement(3330..3918)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     misc_feature    3061..3064
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     terminator      800..843
                     /label="another feature"
                     /ApEinfo_fwdcolor="#85dae9"
''',
        'Metadata1 for part_B', 'Metadata2 for part_B'
    ),
    (
        'part_C',
        '''ORIGIN
        1 tggctaggtc tccaatggtt tcaccggagc tgagatcgca ttaaagtcac actacgaatt
       61 atttgtcccc gcggctccag agacaaattt ctaggagtgg gttcatggag aaccgtgaca
      121 agcatggtag cttcaaatag cgtcatcgga tcttaaagtg attgtgacat atatgggatt
      181 gcgcgtctct cttacacaag ctagtgtgcg cgaggtgctg aacacctatt acgtcaactc
      241 ccgggatccc ggcgggttac agagacataa gcagggtgtc catcattaac ttcgatgttg
      301 agcccgtgcc ctgatcaaac acctctcgta gtagtgggtg ccggccgtgt attattttgg
      361 gggagaggat gcgccctaag caatttcttt ccatcggcca gccattcgac gccggcacgc
      421 aagttgttgc gggcactgtt gcagttcccg tactagagtc tcaagtctag aggtaccagc
      481 atcaggaaga tgggtacttc gtgtggatct tctcctgaac ctgaaacctg aagagctaaa
      541 gcgcgtcgcg ctgccatacg ataagtcgct ctcatttccg gaagacgtcg ttgggtgatg
      601 gcgaacaatc gcgaaacgaa aatcgcggtc ggtatgttga ttcaagacgg gctacgggcg
      661 ggttactaac aggtgcgcag aaatgtgatg ggtcacttga ccatgaactt tcgcgccgtc
      721 gcttctgtgg actggccggc tccctgatag tgcactagcc cgttgatatg cccattacca
      781 tttggcggcg cagtgtgcca gcgggattaa gatgcgtggc gtacaactat ccgcactttg
      841 ctgacgcaca aaggcaactg atggacagcc ttcgggcatt cgactggttg gctgagttcg
      901 cctgcctatt aatatgcttg tggaaattga tttacagcgc ctacttctac ttaggtgtac
      961 ctatttacgg cacagggata aaggctgaca gaaaagtcct ggacaggtca taggtagttg
     1021 tcagacctac aatggtataa cttgagagtc ctgtaccgca catgcggttg tagctagatg
     1081 tccagaagac cgtagctagc gacgagtggg ttgcacgcgc attttgggaa ggggctgcgc
     1141 gagggtgggg tcacgggtta ttccatttct aaatgctaaa ccagcagttg ctttagttct
     1201 gacctagcag gcatgtcata tgtatggcta gattattagt acaatgcagg actaagaaag
     1261 cctgagactt aaatgcgttg gtttttcagc cgacgggccg ggagggcaga tggccgaatg
     1321 cgcactgagg agtgcattgc gagggaatag gcggaccacc aggtcccctc acatactgtg
     1381 gataaaccta cccagggaca actcacgcca tcattttcaa gagtagcact actttccact
     1441 aagtacggtt tgtcgagagc cggtggcacc gtgcgacgag cttgtgggca ctaggcgaac
     1501 cgggttatat aggactgcct tgcggtgcgg gaatacttta ccccgtggag aggccataaa
     1561 agcttacatt tattgtagcc gctagcggtt tcgacgcgct caagtggggc tttatcttag
     1621 agtatcaaat aatgccctct cagcttcatc ttgatttttg tgtttggggt agtgtatccg
     1681 atgggtcacc ttcttgacct tggggagctc catccatcgt accaatactc cgtttctttg
     1741 tgttgcctac cgagtaatcg aataccgcac gagggctaga atctaattta aggtgggtct
     1801 tatccctgtt gtttggacct aacgctgtaa caaatcactt ccgcgacggg tgtcattacg
     1861 tcttctatcc gccaccggtt caacagtgcg tcgaatgcca agtccctatt acccgcattg
     1921 taaagccatt acaatacctc cggtgatata attgattgcc ggttggcaag ggatcatctg
     1981 tacgggccgc gattgattgc tgataaggaa cggcggtaag gaacgtatcg cgtacgaggt
     2041 tgtacactgg aaccgctaaa caagcggtat cgctaaatta ccggcgatta tatcaaggga
     2101 atcttaggta gcacgaagcc ctacgccctg gaaacaggat cggcgcgagg ccgtggacat
     2161 attagcctac gtttttttac acccagattt ccgacgatgt ccacagctaa ggctatccta
     2221 aaccgcaaac acatataacc gcagtataaa gcacgggatg attccttggc catcttgcaa
     2281 gcagtgcaag gcaattctgg tgataggggc aaggcagtga accgagtggt cagcgtattg
     2341 tgaccctgca attgtcatgg cagataattg tttatgaacc tgtttcacac acccatccct
     2401 tggatctcca gcgtagttac taatagtacg atttcctggt gttctgtcac ccgccttgct
     2461 gccttcacgc ctgccccacg ttttgctgtc cactgaatat cgctctaatt aataggtgtc
     2521 acactgtatc aggtaagtcc attccgggtg ctattggatc attcctgaga tgtaccaccg
     2581 atcccatagg cgttaggtct tacatctgga gtgaaaaaga tcagcatcga actataggta
     2641 gataaaggtc atgactgtta acagagaatc gttcctatcg caaccttcac taactggaag
     2701 gcacccccac actttgatta gtaagcccct agcgaatgta gtgatccaag aacatcgagt
     2761 atactttagc tctcggttca ggcatgtatg caatcgacca cccttccgtt cagagtactc
     2821 tcaaccagca gatttctaaa ctttatagcc tccgcccact aaaattgcct tatctcaggt
     2881 tgactcttgg atgtagctgt aatctaacat acacatggag ctgcggcaaa ttatcgtacg
     2941 tccaaaccac caccggatca gtggcgatac tattctcgca aaaaagccat ccagcattgg
     3001 aacaagtgtt aaccagcatt tttctgagtc ccccggaaag aaagttgact acactacatc
     3061 cataaatgct agagtgtgga ttcgccatgt gttacaacgg tcaaccgaaa ggttagcatg
     3121 atacttcgtt gttcaatcgg acgttcgagc gccctccctg tttatcactt gatccccctt
     3181 ttgaaacatg tgcgctgcga cgcagggctg ggagagtcac cggaaaccaa aagcccgacc
     3241 tcggacgtgt agggccctcc gtttgccatc aatggcgatg taggattctc ttttaataca
     3301 acgttcactc aattaacgtc ttggtataag cgtcgaacgt gaccgaggtt gacaccagtg
     3361 tcacggctgc gccggggttt gttgctccgc gtgtacaggg tattggttcg tgagacctag
     3421 cca
''',
        '''LOCUS       part_C                  3423 bp    ds-DNA  linear   UNK 20-NOV-2019
DEFINITION  .
ACCESSION   part_C
VERSION     part_C
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             18..3275
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     CDS             3276..3404
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2486
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b1ff67"
     misc_feature    3411..3423
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    3407..3410
                     /label="another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    14..17
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    1..13
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    3131
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b1ff67"
     misc_feature    2078
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b1ff67"
     misc_feature    14
                     /label="feature"
                     /ApEinfo_fwdcolor="#b1ff67"
''',
        'Metadata1 for part_C', 'Metadata2 for part_C'
    ),
    (
        'part_D',
        '''ORIGIN
        1 accctgcagt ccgtgaatta agaccaatcc actggagctc tatacaggac atcagcgatc
       61 ggtcgagtaa aaccgagttt ctggactccc ttagtctgag aattgtacct tattacgaat
      121 ccggaaatga agtgcgctca taaaatcttg atcgctaact gttttccgct ttttgaacct
      181 acagtcagct atttcgccgc aagcgaggtt tatccccttt tgtgtgtgtc agcgggcgaa
      241 cgtggaggac aattatgacg aatgcctaat aaccaacgta ggcttggctg gtgaggacat
      301 tgcccttccc cttctactaa ccagtgttgt tcagatcatg taccaagtgc agtaatgcta
      361 atcccctgaa ttctgacgtt ggatttggag cgtgattggc aattatccgc tgcaaggcgt
      421 agtatcgcta tctgggaaaa cttagggttg caagcaaggt catcgcccgt ctctagatcc
      481 gacgggagcc tcaccctgca tgaggaagtc ctaatccgct ctaaacaaga gctgaacact
      541 ggatgtctcg cagtagattt gctagaatgc aatgctggct cgtgtgcagc ctcaaggtca
      601 ccatgcttcc gttaaattcc acgctctccg agtctgctgg ttcgggataa atctacgtga
      661 cattcgcgag gtcccggcct gtagatcgtc tcggcaggag agaacacggt tgatccccca
      721 cgcggaaccg atagatgcca ggctagatga tgactgaagg tagtttgtac gagtgacctc
      781 tctagccaag tattttcccg tctcttaagt tatagccgct ctcattccgg gttgtgatat
      841 ccttcatatc cactctctgt aaaatgctgg gttgttcttc tctccacggt cagggaatcg
      901 cctcttttcg gataaacgac attatttcgc gccacagaac ggtttggggg tcgaaggacc
      961 ctagactttg ggtatcccac tcttaccgga tggtaccgct atctccccag ggtccatcgg
     1021 aatggctagc cacgttaccc ttatctgtca gtatcagtct cagacttaag tataccacgg
     1081 tagcgacagc tgtcttttaa tggcccgggc agggagccgg gcccaccgtc catggtccac
     1141 tgtaagggta tctgcaacct tcgccgagct tcttccccaa ggaggtagta ccttaccaaa
     1201 cttccgagtc agtatcgtca aaggggcccc tagggctcac accatcgagt ttccgcggct
     1261 taccagtcta gcctgatatg tttcaggtca ggaaataaga ggtatagccc cgtggacacg
     1321 tactgttcga gcggctagat gtaggttgag ttaagtacag tagacgcgtt ggataccgtc
     1381 gaacattact ccgtctgcca agggtagccg agtacttctt ccggctcggc cattccgact
     1441 aagttagttt cgaattgacg tgccaaagcg tggctcccag tcatttgtcc ttaaattaaa
     1501 cataagtttt ttacctgctc gcgtgccggc cgcttgaggg ggcagcaaga agtcggttca
     1561 aatggggtgt taagaccggt gttgcaagcg gaccatccag agttagcgtc ctcgacgaaa
     1621 cgttaacaat cgtgagatta gaggtcgaat atccccttcg tagggggttt tattgtgcaa
     1681 gatgcgcaaa tgaaccccaa ccttgcgggc ctgtagcgaa caggcgaaaa ggtccaatac
     1741 ggcgcccgat gcaccgtaaa acaggtcctc taactgtgtg ttctacgctc tccggatcct
     1801 gtatataaga acgaatcccc ttttcctagg gcccggccgc gtagacccag tacacttgac
     1861 tttcacgaag atccatccta ctcccatacg cttgagagtg ccaacgtagt ttgtaaccga
     1921 ccttgcctcg gctagaaagg cattttgtcg catgtggcca gcctgtcagg gcgtgctcta
     1981 ggcttgacga ttagtgctac gggatgattc gtaaccgagc tgacgggacg cctctggaat
     2041 tcggtctctt tcgatggtgc ctcctctacg tagatctact ggataaccgt ccccaatatc
     2101 cgcttcccaa cgtagattgc cgacaggcat aagcttcggg ggcgcaaagg ccgacgtccg
     2161 cattgcagtg tagctttgtg agcaggaagt gtgatagtct ttcgattatt aaagtctgag
     2221 ctgaatgaaa aaaggtccaa cgaatggagg acgcgaagac atgggtgtct cttatggccc
     2281 gagcgggagt aatggcggtt cgtacataaa ggctgaaagg attctggcgt tagctgtctt
     2341 acgttggatt ggcccttcaa attatcgatg ttagctgatt cggtgtaccg ggcgaggaaa
     2401 gcgctctcag aacaacttca tatacgaggt tcgactataa tggtctaagc tcctgggcta
     2461 gtctcaagaa gcgggtacct ttagtagcac gtatcgacgg caaagcaaag aataaaaact
     2521 tggctttgca tcgtgcaaag atttctaact aggttgttta agggctggta tctatgtccc
     2581 gctataacag cgcgcctaca gtagaagttt aaccatgaca tacctttgaa gtgttcgtat
     2641 cacacacaag gaaggagcat gtggacacca ctgagctttg ag
''',
        '''LOCUS       part_D                  2682 bp    ds-DNA  circular UNK 30-OCT-2019
DEFINITION  .
ACCESSION   part_D
VERSION     part_D
KEYWORDS    "Source:Subcloned from vector Andrew" "Sequence" "Freezer" "Shelf"
            "Box:Mobius box" "Concentration" "Use:Mobius".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     rep_origin      complement(262..850)
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2642..2671
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b1ff67"
     terminator      1991..2034
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     CDS             2054..2566
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     misc_feature    1919..1938
                     /label="another feature"
                     /ApEinfo_fwdcolor="#f58a5e"
     misc_feature    2050..2053
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    2576..2641
                     /label="feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     terminator      10..67
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    143..162
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#75c6a9"
     misc_feature    2675..2678
                     /label="feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     CDS             2618..2641
                     /label="another feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     CDS             complement(1148..1807)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     terminator      1033..1127
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
''',
        'Metadata1 for part_D', 'Metadata2 for part_D'
    ),
    (
        'part_E',
        '''ORIGIN
        1 tggctaggtc tccgctttga tcagcacgcg tctcagagtt tcagggggac ccaaattact
       61 ggtcctcaat tgggacgcga ccgcatctcc ccacgaaagc ttatggggat tgcccactgc
      121 ccagcttcaa atctgaaggt tcggcttatt gacagggtct aacacgcagc tcaactgctc
      181 gaggttagag gcgtaatgac gggccccgat agccttctac gattacgtcc agcaggaacc
      241 ccacagcttt tccctacgtc taacaccgtg aaagcaaaac tgtctgccct tttacatggt
      301 ctttaaggaa tctctgcact tatgttatta ggtatgagac ctagcca
''',
        '''LOCUS       part_E                   347 bp    ds-DNA  linear   UNK 20-OCT-2019
DEFINITION  .
ACCESSION   part_E
VERSION     part_E
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     3'UTR           21..330
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    1..13
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    18..20
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b1ff67"
     misc_feature    14..17
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    331..334
                     /label="another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    14..17
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_feature    335..347
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
''',
        'Metadata1 for part_E', 'Metadata2 for part_E'
    ),
    (
        'part_F',
        '''ORIGIN
        1 accctgcagt ccgctcacgg accgcaagga cgggctaatt aggaggcaac gccgatgggg
       61 ccgcagttca gcgctgcaat gtttgctgaa cagggatgtc acgcatactc gtctacaccg
      121 cccgtgaccc gttatcacgg ttgaagtgtc gaggactagt gctgccgctg cgtgagggac
      181 acaagctgct atttgtccat acgccatgtg ctccgagctc atgctgccat gagacaatga
      241 gacatgtcgc caataatcga gtgacgagtc agaatgacct ggctccgcat aaccgttcaa
      301 agttattgac aacgcatctt tcgtagttcg tgcagcagcg gtctttcttc tatagccgac
      361 tagatgttaa gggactcctg gataccgcta gttttaccct ctccaggaag ccagcgaggg
      421 cgtgccgcaa gtcccaatag ataccgggca tgatcaaggg gccctgtgct ctgagtctgg
      481 aggcgacagt gcgctgcagc tcagaggtgg ttattgcgaa ccggcaccgc tggacagcac
      541 ccacggggac acgtaagtaa tttagggtct gggccaacgg ctcagcgcca gtaggattaa
      601 caaactcgac taatcaatgt gccagctact tccgccgggt ctgacggggc ggcacccatt
      661 accatgtgta ctgaataggg attccgagcg accgtaaagg cgttcctaag tgtcatatac
      721 tggcaactag aggcactcac atcggggtta agggcccacc gtaatggcca cgcaggatac
      781 caattggccg ggtgaggatc tattcacgcc gatggggagt tctaagcccg agttattggt
      841 gctagtggct tggacctgtg tgtcgatgat gcgcgatata gaggcgcggg actagctggg
      901 gtccacaggt gctatgttgg gtgcgctcta tggacgtccg gaaagagact aatgcaacgg
      961 tatggcatca agcgcgaccg atggggagac tggaaattgt gaaatagtgt actggcgcga
     1021 tcattaaata ttgcatacgc tgtccgttat gacctagagg ggattatttg aacgagaggt
     1081 cttgggaact gcattggaga tggtcagatg gaatgcgaaa tgtatcacca cgggcggccg
     1141 aaggggagac agactatccc tagttggctt cggtaaccgt aatcagatgg aacggcgggt
     1201 aatgctgata ccgagcgctc gtacggcgcg tggtgatgac acggtctgat gtcgtagcaa
     1261 acgagccggt cctcgtacaa caggctcaca tctcaatcac gcataccaac cgattacata
     1321 acgcgatcta tatttgggga actctactta ctaccctgtc tgcagagtgc gttctatgat
     1381 cctcccctac gtgacggcca acgattagtc ggcctaggtc taacggataa aaggactccc
     1441 cagtacatgt gaccatatga gccaggccgt cgagcgcgac cactccccag cctatagtga
     1501 ggaggagcgg tattcggtta aagtttagct agcgactttt ttgtcaccga agtagggacg
     1561 ggcatattgt ttacccttaa agcggggatc aatccattgg ggtcggcagt cataaaagga
     1621 tctaagcccc caaagcgcat ggtaaggtac tcccgggttt cccgaaatct aggcaagtct
     1681 cggtgggtgt cagcctgacg agtagggcac gttaccgagg ccacggtgct tgatgactac
     1741 gggtgaatcg aaccgatcaa acgcacagca tctaaaccct ggtgtactgg tccacagggg
     1801 gacaagtctt tgaagtgtcg cttcagatga cggccgcggg ccaatctgaa ttgttagaca
     1861 ccgacagtag ggtgtgcaac tcgctcgggc aagagtgtaa aggcacatcc tccccgaggt
     1921 aagttatacc gctctaacgg cgcgggcagc ttttcaactc aacacttccg cggttcagtc
     1981 ctgagcatta ggctgcatgt ctatcacaag aggtgcgggg aacgacaagg gctctggaat
     2041 tcggtctctg gtaagctccc ataagagcac cacttcgtgt accttgctaa ctcccttcat
     2101 ccccacccgc gaaaacttaa gacgtcaccc ctgttattat cccgtgcact cctaaagtgc
     2161 gatgggcaag aacgcaaaca attgtctgag atttatatgg gcggcggaca tagtcagaga
     2221 gccttatata atctcccacg tcttgccagt ccgagttatg agaaacccgg gaccgatgac
     2281 gctattacgg tcgcccagcg aggttcaaca gagcgaggca taggagtcat tacgcgtgtt
     2341 agcttcaaca tgcgtggaag tgaaacggat ctataaaacg gcgtggagtc atcagcttga
     2401 gcagatctaa cttactcgcc acgcgcgcaa atcgtctttc gctgtccacg ttatacgtaa
     2461 ctcgcttgag
''',
        '''LOCUS       part_F                  2470 bp    ds-DNA  circular UNK 20-OCT-2019
DEFINITION  .
ACCESSION   part_F
VERSION     part_F
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     terminator      1991..2034
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     terminator      10..67
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     rep_origin      complement(262..850)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     CDS             complement(1148..1807)
                     /label="feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     misc_feature    2054..2462
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
     primer_bind     complement(2351..2371)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#85dae9"
     terminator      1033..1127
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
''',
        'Metadata1 for part_F', 'Metadata2 for part_F'
    ),
    (
        'part_G',
        '''ORIGIN
        1 tggctaggtc tccggaggaa cttctatggg acgaaagatg cgcagctcga caaatctcaa
       61 attggaaacc gtctaaggtc tttgagccca cgctaggcca ccgaatagtt gcgagctcgt
      121 gggcaattgc tggctagcgg tttgtatcgt atacgtaaca agagattgcg ccggattcat
      181 tcgtgtgggg tgctttcaca tcggaatact caagagtggg ggtttgcgct ttaatatgac
      241 ggcctcatgc accccaaaaa taagagcgcc aatttcatcg accacattcc ggacaattct
      301 gacttccctg tcggatcgac cctcgctcac aactcctaga actccaagac ggtaggacgt
      361 cagccgaaaa ggaggagtga cgacgggccc tcctagctga agcaatgggg ggcgtcgacc
      421 ctagcctggt gtgggacttc ggagggtcgg gcgtgtggcc attcctgcac gaggggccct
      481 tgcagtatcg aactttgacg agaaagtggg gggcgaccta taaacatagc ggagttcaat
      541 actcctgcaa gtgcatgatg acgttcagcg ttggcacatc gaacgcgtcg ctacacttgc
      601 tcaccggaag gcaaaataga ccggcagctg gcccgcacgg acgtctagac tcctacgctg
      661 tcaaatgcac gatgactatt agcatgcgga attcgaggcg gccggtacat cgagacacgc
      721 tggtcttaat acctgtgtta tgtcaaaaca agtgtctcgg cggctttgta ctacgctgtg
      781 aatgcgcatg atgttgcgag aagctgaacg ttggtagcag tctacaacag aaccgacgag
      841 ctacggacgg gttacagcag gatcccctct taaggcatta ttcgacgact ccaggttcta
      901 accgatcagt aattgcctcg gatggtcgta cgtgttaacc gagacagcaa ggcaccacat
      961 agacaactgc atagggcgcc tacaggtcca caatcagggt ggcccaagat cctcaaagac
     1021 ttgttgcagt ctctgcttat tccgtttaac acgttgtgtg agctctagct tactacctca
     1081 tcggtgtcgg gagtgattta ggaatgtgag acctagcca
''',
        '''LOCUS       part_G                  1119 bp    ds-DNA  linear   UNK 25-OCT-2019
DEFINITION  .
ACCESSION   part_G
VERSION     part_G
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     GoldenGate      1..13
                     /label="feature"
                     /ApEinfo_fwdcolor="#faac61"
     GoldenGate      1107..1119
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
''',
        'Metadata1 for part_G', 'Metadata2 for part_G'
    ),
    (
        'part_H',
        '''ORIGIN
        1 aacatggatc cgttacacaa ggtctactcc gcgatgtggg ctcaccccct agatcccttt
       61 tcctaccgcc cgatcgcaaa ctacgcgtag cttcacggct ctctgtctgt cggctgacct
      121 gcggctcatt tttcatccat ttagtatagg ttgcaagggg tgtgactgat ttctctacat
      181 aaagaagcct atagtaccac atgcaatgcc gatccggtga gtgctgtagt ccttattccg
      241 gcttggtggc ctcctgccac ttaggatcgc aacggaattc tagttcctaa cgtcccttct
      301 acctgattaa gtgagagaaa acggaaccaa cctaagacat ctatgcatcg atgttttaca
      361 atgacgagag tgctccatgg ttcatactgg ggacatagtt tcgggcattt aggcccgcaa
      421 gagttcgggt atctagtttt tccgtacgaa aaacaccacc ttagaacagg tcacacgagt
      481 acagagagta aacccgtaat gttttcaccc attccccgat ttgactgcaa acaaacaagt
      541 ctccccgtat cgcttcattg ttttatcaac gggggcggaa tcgtagacca ataattaacg
      601 acaactctag gcttttcaac cgacgggggc aaagtctgaa atgcctgaga aaggcacaac
      661 gacggtggga ggggcccctg caaatattgt tctagcactc caggattcac accaacctcg
      721 gggccagcca ttccggtact tttggttata cggcaggata ccacgtccct aagcggagta
      781 tacacggagt tgccttagcc aggtacccgc ttaaccgtcc aactcccctg agtggtacgt
      841 ccatcacgcc gcgtaaacgg tttccaaaag agtaagctta ggcattcctc acctgtcacg
      901 tccctcttaa gatccctgcc tcgaactggg atcttgtaaa aatgtggcca tgggagcagg
      961 agttcatgtg gggcgacggc gatcggacgc ccctttaaac cagagccctt gacgctagga
     1021 gatcagtaga gttacccgga atgagatctc taatctaagt ttgatacgac gaaagcggcc
     1081 tgacgggcgc tactcatgct cttagaatcc gtcacttatt ctgccacaga tccgaggtac
     1141 tgtggatgtt atttgcggag actggcctcg gattaggggt actaagtccg agatgtcatg
     1201 agttaagccc aattcactag ttgcattgtc aacgagtgga cctccaaaac gatgttaggg
     1261 tcactacccc aagcgagcac ctccgctagg accacacacc cccactttca aggtattttg
     1321 ttcgcatcac ataccgtctt ccctttgcgg tctaatagtg aagcttggga caaggaacgg
     1381 ccagtcgcag caatactacg agtgagagtc ctaggcgcaa ctacgtaggt tccgcaaacg
     1441 tgttgtggat tactgttgag gtttgagttt atccggaggg ctgtgtaaga attaccagct
     1501 agtcaaaagc cctgcatggc ctgatctcat tagatacctc ggccgggccg aggcaaatcc
     1561 gtgaacaaaa accatgcatt tattctatca tagaaacatt tgttacacat ctactggccg
     1621 tgtgccgtga acaagaactc aatttagtta tcaagggact gctgtaaacg gagccgcgtc
     1681 acccgcgtgc acacgtgtag tgcttacgcc ggcccccgtc cagcgacacc gtagtcaagt
     1741 aaaaaatgcg tattcacgac ctcacgtacc cgtttcggag ggtgccttga cgcatagagt
     1801 tgtctgtgtg attggaattt gtaaggggtc cgcccatgta aaatagagcc catgtcttac
     1861 gaggacatga ggaaacatag gttctggtgg ccttccaaga agctgcctct acactccttc
     1921 tctttaatca ccgaaaacct taccttgaga ggacgtcatt cgtcaaaaca aaatattggt
     1981 ggcaaattaa gattaccact ataggcggtc tcaaatggaa gacgcgttta gtcttcggcc
     2041 cccagacggc accggcagag ccactgccaa tgccttttac gctgaactct tcgactttgg
     2101 catgtgcgac tgactatacg ccgatggccc taccctgcag cggaaaggca cagggccaaa
     2161 aacagtggaa gcgcaccgct agagctttac taaatccttt agggtgaacc ggtgttcagt
     2221 caacagaaac atatgactag tccttaatgg cactactaaa ggcggcccat ctcttttgta
     2281 ggtcacgtct catcgtttag caagcgtccg ccatccgagg agtatcctag tacgtgagag
     2341 ggtgttttac ctcatcacat ggtttttggc tagatgtggc agatgccagc tgatacctta
     2401 caacccgcta taacggtttc tggatcgacg atacaacctg cttcggcata tcataatggc
     2461 ggcctcatag gacttccaat catggttatg cgagtttgat gtgttgaaga ttgtgagggg
     2521 aaaccttttt ccgagtacgt tgaatctgcg actgatggtg aacaactcat cgagagggat
     2581 gtgcgggcag tctccttgat gcgcgaaaag gagtcgcgac cgcgtgttcg tgacgcactt
     2641 gatgctagac cgtaggggct cattatctcc tgaaacagta gctagctttt ggtgtccgag
     2701 tcagcttagg cacagtaatt gaaaagaagg aaatgcgtga aaactaaacg tagacactcc
     2761 ttgctaagaa ttccctggtc tttttcgagc atagactacc gtaggtaaaa ggcttgtcgt
     2821 gtcacaaata gacttgatct gttcggcata tacacatcct actgcttcag taaccaagga
     2881 cgattctccg gccattggac tgtgggacca cccgcgttca aacggaatcc tctacgaacg
     2941 gttaagaaag gagcccactg cttgcggaga catatcgttg ctgagtggta caccggtgag
     3001 tcaaaaagta ctgacactac tggtggactc tcgggtacga tgtgtggggt attggcatac
     3061 tatggcgcta agcgtggaat agtgtttcca atagcccttt acctaagaag cattctgaca
     3121 catgggttgg acttttacac agacctataa ctcgtatccg gcctggtgat gtatgtccgc
     3181 ttatcactga atgcaggatt ccgccaattc cgctcgcggg ttctcggaca acgggccatg
     3241 gagatggatc aatgagccat gttgattacc atttcactac acacggttag tcaattccgg
     3301 agatcataat tcgacaattg aatatgcgtg atacgctttc tgttaccagg atcgttgcta
     3361 attagtttcc aatgtcacag tgttgggcgg cgaggtatca cggcgctctg ggtagaatgg
     3421 cggtgcctga ctggacctca taactggtaa catacctatg accaaaggag gtcagacacc
     3481 gatgagtcgg tgacatgggc ggtcacagcg tgccctgcat agtaatacag ctcgtttcta
     3541 gacatttgtc ctaattcgcc acacaagtac gagagcgggg gtaacccatc accatcgtgc
     3601 ccgtaggtac ctaatcttca tcgagcgcta ttgttaaagt ggactcaggg atcggcatct
     3661 cgtgtagctt tgagaccggg gcccgaggat gtatgccttg catggaactt gactaatctc
     3721 agtaatactg gccgtgtatt ggtgtgtctt ccttggcctc acatgggaac attaccacat
     3781 ataacatcat ggggttcagt ctattactaa cgaagtcact aaaacgtact cgacggcgca
     3841 tgcgagtatc gacgtagcgg ttggttacct taagctgact tggtctatgc atacactgat
     3901 tagcttgtag cggtggcctt cgacgtctgt atcgtaggcg ttgtttgttg caaacacaaa
     3961 cgacctagcg tcaagacccc ttctagagaa agagagtgac gactcttgtc agttcaactc
     4021 tgtggtacgc cgctgcgaga tgcgcgttgc cacggcagac cccatgaact cgttcaagct
     4081 aaactctagt tcacgctaca tgtgtttgtt tctgcgtacg caaaagtact aggtacagct
     4141 atcttttggg ccccaagtgg gaccacactt tcggtcttcc tagcaaggct ataaaaccgt
     4201 atggcatgag tccttgatct gcatttcgtt tagacagagc gagatg
''',
        '''LOCUS       part_H                  4246 bp    ds-DNA  circular UNK 07-NOV-2019
DEFINITION  .
ACCESSION   part_H
VERSION     part_H
KEYWORDS    "creator:militzis" "marker:SmR".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     primer_bind     1966..1982
                     /label="feature"
                     /ApEinfo_fwdcolor="#a020f0"
     terminator      3774..3860
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     rep_origin      complement(63..651)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     CDS             2015..3667
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     promoter        1987..2005
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     CDS             complement(744..1535)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     misc_feature    3667..3670
                     /label="feature"
                     /ApEinfo_fwdcolor="#ff0000"
     terminator      3952..3979
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     primer_bind     complement(3694..3710)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#a020f0"
''',
        'Metadata1 for part_H', 'Metadata2 for part_H'
    ),
    (
        'part_I',
        '''ORIGIN
        1 accctgcagt ccgcatccaa gcttgctggg gcctctatcc gggacctgct tccgtacccc
       61 gtccgcctat caggctaacg actctgcctg ccattcaggc tcaccggcaa gttggagtga
      121 cctcttcatc gatacaaata gcgcatacgc cagacatacc gcccacgcac tcgaccataa
      181 aattccccgt cacggcgcgc acctccacac tccgtggcct aaggaacgtg gcagccggct
      241 gggccacttt ttagactacc atcggctcat ttttagaagg ccaccggctg tccgtcaatc
      301 gtacccgttt aggcctcctt aacacagttc cgaatagtta cctcaagatg cggatacagg
      361 atggctcccc caggggtcta gttcagtctt ttgtcgcctg cgcttgaacc tttatgcagc
      421 gatccacaaa cggaacaccg agatttgtat ttgggagggc aggctttgga ggaacggtgt
      481 atgtactttt gcttggtacc ctaggtaccc cacgatctaa tcggtcttcc tgacaggtga
      541 cagaggaggc cgtaactgcg accacaggct aactgcggcc ttaacttatg gtcgcgaact
      601 gcataagggc gtcgccgccc gttgtcagtg cgaattctgt tagtcgtcgt gtacaccttc
      661 cggtctaaca tctcgtagtt actagaagag cggatgggtc ctgtaaagag agcaggttcc
      721 cgtcttctct cactttaccg gtactggtgt gaccgggtgt gagtgatact gtgattccgc
      781 tactattgct atgcctgccg tcgtagctag atacgatccc agaacatttt gggctgattt
      841 caagtctctc cgggattccg taacgacggt tctcgaactg aatagctttg agatatcgca
      901 atattccttg gttgactctt gtccccggga cacaagtcgt gtgcaatatg tatacgtcag
      961 tatgctcgag actcctaagg ccgcatatct atgctatcac tattgcctat gcagtaactg
     1021 gctaagcttg tggggtactt gcacatgact catgtcaagg tcggaggatt cccagacagt
     1081 tggcatcatc agtgcgttca aggcggggag cgaccagcca gatattgatc gacgaggtgc
     1141 cacgaagtgt tccaaggtta tttttagtat gttacatcca tcagcgggtc taggccatac
     1201 cggtctatat tatggtgagt cgtataggct gtaatgccgg ctgcacaatc tacggccgat
     1261 tgcacataaa atggtcaact aacgacaagc tctgattctc aacttggata aaccttatca
     1321 caaaaggtcc gacacggttc agccgaggta aatcaaagtg attatgactc caggagcgac
     1381 acaaagttgc actgtctact gttggtcatc accgtactgc aaacgggtaa agatacttag
     1441 ctttgttagc tattgcagaa cctaattgct ttcctgccct aaggcgatcc ggatcctcta
     1501 tttatcaagt taattatcag agcttagctc aagtcaaagg tcttagtcag gtatggttta
     1561 gtcggcttat ctcgctccga taaacccctc gcgcctcgga ttcattcacg cgtttatatg
     1621 gtgaggggcg acccgtaggc tagtcgcact cctcggcaat taatcctaca gaaacctaga
     1681 catgggagtt gcgatcctcg agatacggcc tgagagggtc ggcaaagtgg gtgcttcatt
     1741 ctcttcggcg gagttcgcag gctctagaat gcctggtttt catctcgaca taaattacga
     1801 tgttatctcg tgcattatta ccccttttct gcaggtttta agtcgagtct agtcttacta
     1861 acgctgttct tcccgagaag gtgtccaggt tagatctcgg catattttca ctggcgcgct
     1921 atccgacagg aagacagcac ccgttaagcg ccccagtcac atcaagatca gtcctatgtt
     1981 gactcttgct ggaggtcttc tcgaatccag ctaatgttgg attacccgcc cctctggaat
     2041 tcggtctctg cttttattat ccggcatgca tgcggcgttg catcgcttca tcgtcatcga
     2101 catatcgtta gcgggactta tcagggtccc atcacggctg tgacagctgt tcataagtgt
     2161 aatagcagtc accccaccaa atgcaagact gttcccaaac taagtgctag ggaggttggc
     2221 tacctatcgg cgtggcagtt gaccgggttt ttcttagtac ctgccttcca atctccgata
     2281 aaagcaagat gacttcctgc ttgcgctaaa gtctgaattc gctaaggggg cagattcatt
     2341 accagagact gccaaagtca accgccgaac gatactggaa tatatggagg tatgag
''',
        '''LOCUS       part_I                  2396 bp    ds-DNA  circular UNK 29-OCT-2019
DEFINITION  .
ACCESSION   part_I
VERSION     part_I
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     terminator      1991..2034
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     terminator      10..67
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     rep_origin      complement(262..850)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2054..2388
                     /label="feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     CDS             complement(1148..1807)
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     terminator      1033..1127
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
''',
        'Metadata1 for part_I', 'Metadata2 for part_I'
    ),
    (
        'part_J',
        '''ORIGIN
        1 ctcatgggag tgtagcgcta gccatagttg acctgtcagc cggtcaaggt tggtccgaac
       61 ctagccagtt agacgaccag ctcccccact attgcacagt gcctgtataa gcacgtcagt
      121 gtctggcacc actgtgcgcg gtgtggtcct cgggacgatc ttcctacggc tacccatgcg
      181 ctagtgatta cagacgatcg ttttttcttt ttctgtggtg caaatacaca agatacgtca
      241 agagtcctag cagcctcatc ttgccaattg caggggtacc gattcaccca tggtatcagt
      301 cacgcaaaaa aagacgacat cggagttgct ggcgactggc aaaagaaaac atctatctgc
      361 cgtggcgttc gtctagcgca gcgagccgac tgactactcc gttatctcaa atttagtttc
      421 gtaagatctt cacgggagcc ggccagctga agaatctata ggtctaccgc tagaggaagc
      481 ccgtttacaa ccgacgtact agtccgacag caatcgtcag catgccaatt tattgtatca
      541 ggttgagaat gcaccactac ggaatcgacg tgccactggc cgcacgatac tagtgcgctg
      601 gtaaacaggc caggagctct aaatctggag tggtgtgaat ttaaccgcgg aacgattgtt
      661 gcacgcggtc cctgtgactc tgttgggagc actgtccccg gctccggata agtgatgtga
      721 acgaaaaggt gggctgacca gagatagggt taaggacggg cgatcgtagt ggacgattac
      781 tgctgatctt ccgagcactc tgagtcgcgc ggcgtatcgg gagttacccc gctcgctaca
      841 catctgaatt ccgccggata ctatgggtgc agagtggatc accttttcgt gacccgtata
      901 ttcgtcgttc tcatcactgt cctccgcgat ttatttcaat cgtcggtttc gccagatgca
      961 ttcctgggga agtcatcgcg acagtgcgct gcgcgcacag ctcttcttgg tttattcgag
     1021 attgtggttt atcaacatgg tctatctacg gactgattag aatattccca cgttaagccg
     1081 acgtgcattg gtactggctc cttaattaag cgcagagggt aggcgccttt tcttgtgcaa
     1141 accattatat gaggaactag caaacactat ggcgcaatat gtccacgtga ctcgagcatt
     1201 ctgtgttgcg gtagcgattt cactagtcac aaacggatag attctacggc gcacccggtg
     1261 catcgttgta tgttcagcaa tttggtcata tcccgctcct ctattcaatc tcttttttcg
     1321 tcaggggggg cacatacttg gagcagctgc gcctaacatt acaataagca gtcggagctc
     1381 agtacaaaac ccttatctta gcggtccgtt ctggaattcg gtctcaggag atatgcaggt
     1441 gtttacggct aggcacttac gcgaaatttc agctcgcggc aaactgtgct tatgcgactt
     1501 cagctagtcg tcatcgagtg ctgacgcgct tgttctacct taaggtgccc gctatcgcat
     1561 tgtccagggt tcggctatga gggtcggaac cagatcgcac ctacatgatg tcacgatact
     1621 attattcacg tcgttagaga acggacggtc cttcgtcctc acaggaaaat tattagaaat
     1681 ggtcaaccgg gtccacaaac gggcaattgc cactttggta gtacggttgc acataagatg
     1741 atatttacgc cttctgcgtt caggttacga tagtcgttga tggggctcgc tccgggtgaa
     1801 gttgatccag taaaatgggc gtcgatactt ctcctgtgag ctaagcgtct aaattcattc
     1861 ctagtctcgt gacgagaaca gaaaagacta aaccacgcgt gatggacctt ggtaagcgct
     1921 tgtcggtcta acttgtctca tctcacacgg tatgttcttg tgagcgctga ccagcccgtg
     1981 tccctgacgg ggaattcatt ctaacttacg gcttgggcgt cgtgatgttt cggaatccag
     2041 gatattacac agtagatgta gcttgtggcc gaacgagtgt atttcctgcg tcggattacc
     2101 tcatttcttg ttaatctgcc tacaattaag ccctttcagc agccttccct ctcgtttata
     2161 tttcgacagt caccgtgaca gtgtgatcga gctgtgggaa aatcacaata tagacgttag
     2221 ttgcttcgac agaaacggcg aatgtacgcg gttgtgatgg atgagacagg ccattcaatt
     2281 ctattcagcg aagttcagat tggttgctac tgctaaagga tgcgcatttg tcctttctgc
     2341 gcacctgcat atacccctaa gatgccagcg caatacggca tcgcagggca agcaacaccg
     2401 cacagcatgt ccctaaaggg aagatttaat acggactcag tcaagttgtg agaaacgaca
     2461 atacattgag tggaatgccc tttggtacgt tggcgataga tctaacgaga aacttcaacg
     2521 agatctagag tatcggactg agatgcatca tggctaggga ttgctgaatg gatatctttg
     2581 ccaaaatgag ctcgactagg aatctcagac tgccagctac gaaaagctgg gggatattca
     2641 tttccgctct tgtaaggcgc acgatgaccc tgtggagcaa aattcgacgc ggcaacgtag
     2701 tcaaagtacg gcctatctgg gacggagatg ggattgtttc tgtattcatg ccggcctctt
     2761 gagattgtaa ggatacccgg ccacagtctc gtaagaaacc cgatatcgct caaacgaatt
     2821 gacatcagtt agggcttgcg ctctagcaac cacgcgtcgg gtgttcgact gaatttaaga
     2881 tccagcgcaa accattcctt aaaggcggtc caggtgaaac cggctctgat aatgtacagt
     2941 ggattgtcca gttcgcttca aagtgcaatc aggcaaggta ctggcaagtc accgcttggg
     3001 gccggtttca gtggagcttt atcgacaaaa ggccgggccc tttggcctag agctcacttt
     3061 tggtgtcgtt cacgagtgag aactttggct cgcagcttga atgcttgtca cgtttagggt
     3121 tatacaacgt ttctataaac tacctatgac aaatcccgat cttctgaaac tccataagga
     3181 ctaagagagt ccgcgacagt ttcacttggt cggtgctatt gactatgtca cgcatgcccg
     3241 gtaaagcttt taagcaaagc ctgatatatc acgggaggtt actcaaagca ttttctaaag
     3301 ggacagcgga ctcacattag cctaatataa acccgatcag gccatgcacc tcgattcacc
     3361 ggaaggtagt accatgcacc ggatatttgg gacaacagtg aatagtacat cgtaaaatgt
     3421 caaggcctag ctatttttca ctggtgtctt atctctatta caacacacac ggatgtcccg
     3481 agcctacact ttcaaatata ctcaaggttc gctgttaagg gtaatacaag agtgctcggt
     3541 ttagtatcca tatgttggta caagtgactg tgctagccgt aaattaattc ggcttccatt
     3601 ccagcttggg tgtttagggt ctagacggtt gaaaccagaa agagtacaga caaaaccgta
     3661 gctctcccaa ggttgatcct ccagacacct acccacacta gtagcgcagg cctagaggag
     3721 acgttaggag cggaattgta cattcagtat cgcattaagc acaaagacag acaaacctag
     3781 aagccactta tcccctcgta ggagtacaac ggcagcgctt tcgggaggag ttgcctaact
     3841 acgctctgga cagcaagtcc cagaatgaag aaaatgactt ggcgagaacc caccactaca
     3901 tgcagggtct tgggaccggt cagcacagct tatccatgcg gcaattccgg acgaccgcga
     3961 atgaggatca cgacactgcc gatgttgcac cattgacctt ctcctatgca acaagggctt
     4021 aagatggtag atggtaatag agtagcctgc gaccttatgc ggtaatacac aaaaaccggt
     4081 aaacagtgcg ttgatgctgt gtgctggcgc atttgtgtac cttgctgggc tcgtataaaa
     4141 ttctgtcagc aggacgttca cctcccttct gctagtcgct aggcacaaac ggaaaaggtc
     4201 taatcgcgcg cgtgctcgtt agagtactgg ggcctcccac cacggatgat cacctcgccg
     4261 atacatttgg ccttatgtct caataagacc tcatggtcta tcctacacca tagatttgga
     4321 aagtctagtt gctgcctagc aggctgatag caagcgtcgt acgaggatcg aggaagggca
''',
        '''LOCUS       part_J                  4380 bp    ds-DNA  circular UNK 31-OCT-2019
DEFINITION  .
ACCESSION   part_J
VERSION     part_J
KEYWORDS    "Source" "Sequence" "Freezer" "Shelf" "Box" "Concentration" "Use"
            "creator:SynthSys Center" "marker:SmR".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     rep_origin      join(4298..4380,1..166)
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     promoter        1442..1476
                     /label="another feature"
                     /ApEinfo_fwdcolor="#85dae9"
     misc_feature    2345..2362
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             complement(2550..3341)
                     /label="feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     misc_feature    1187..1210
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    1366..1396
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             complement(2550..3455)
                     /label="feature"
                     /ApEinfo_fwdcolor="#993366"
     misc_feature    2421..2449
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -35_signal      3766..3771
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             167..1186
                     /label="feature"
                     /ApEinfo_fwdcolor="#993366"
     misc_feature    2526..2549
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     gene            167..1186
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -10_signal      3744..3749
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             1503..2222
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    1409..1426
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     -35_signal      3589..3594
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     terminator      2227..2298
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    1319..1343
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    1370..1394
                     /label="feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     misc_RNA        3625..4177
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     terminator      2314..2341
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    2364..2474
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_RNA        complement(3628..3735)
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    4264..4272
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     RBS             1485..1496
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2353..2356
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_feature    complement(2476..2516)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    3539..3562
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    1427..1430
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     -10_signal      3610..3615
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     misc_feature    1211..1258
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     gene            complement(2550..3455)
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     rep_origin      3589..4177
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    1414..1419
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     misc_feature    2388..2411
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b4abac"
     terminator      1409..1411
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    2357..2362
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
''',
        'Metadata1 for part_J', 'Metadata2 for part_J'
    ),
    (
        'part_K',
        '''ORIGIN
        1 accagacagc ttccctgctg cttagagtgc catggggaat tgatctgtgt tcaatgactt
       61 taccataagc gcgcatcgtc attgcatcac gtgctttatc tctcgtagct agtagggaac
      121 aagagtttgc gaccactggt ggaatttcta gccgctatgt gaaaagtcac tagtagtatt
      181 aatagtcggc acatcgcgta cgcagtgctg agattcccct cactgttact tgagctagtc
      241 gcgccttgac agaccgcctt tcgaaggtgg cagagtgcct aatatattgc cacgttgagc
      301 gtactccttc ccggagtttc atcttaccgc cgggcgcgcg aagctacgtt tttcgaattt
      361 actaacgcac gctgacggtg gaagggcgat tggactggtt agattcaggt atgccttaag
      421 cgtactttac aaccagaggc tttatttgat tgtctaatct tcgtcaaaga acgaatacca
      481 gcgcttgagg ggtgactaga ctaataggag acaccgtctg tatcactcgt cggctatggc
      541 cgtactgtga acgcgacggg ccctaccagg ccacagattc ctatgtacgt tgcaagggtg
      601 ggctttgacc gttgcaacgg cgcactaggt ggcattcttc tctgctgcca atgatccgcc
      661 tcgccgagct aggggcgaag caggttatga cgcagcggaa cggtagtact atctaatatc
      721 aactaaaatt gtttctcaaa ggttcaaact agtattcttt attagaaacg atgatggcat
      781 ccggaacagg gttatggccc aggtgctgga tgaagccttg cctcgaggga cttagtcgcg
      841 tcgcccctct accaaggcgc taacggctag cgagatggtt gaggcgggag ccccgccact
      901 ttctagttga attcagagct acgtatcctg atgctaatcc ctcgaagaac tttccgatgc
      961 agctactatc tctttcgtaa agttagaata ggaacccggt aacggtggaa cgtcccgtca
     1021 taggaggttg tctcctacac tatctaagta cgtcagaaga tcgatttcct cagccgctga
     1081 ccccgtccgg catgaccagc gccgaacgta taggccgtat tgctgccgga acgactgctt
     1141 caaaggtaac catagccgtg gtacagtgac aacccactta gctattaaat catgcgcact
     1201 ctagacctct tgcacatacg gagcaactat gtatctaatt ggcacgagac atagaggagg
     1261 gccaatcaca gtagtttact ggtcgaattc cacctgcata tggagagaca ctccttcatc
     1321 aacagcagag gttacatacg tttaacctga ccgatggttg cccatcgaaa gcaaggaccc
     1381 acgggtaggc agctgttgtt tgatatgctt ttgtccccgg aaggcgccga taaaattata
     1441 cgccgtcccc cgcgtccgtt acgagcttcg tgatgttgcg ggctgaagtt ccctacgagt
     1501 caataggcca actcgaatgg ctaaatggtt gcaacccaca ggctcaactt tgaaggtaat
     1561 acgccaaaat ccaatatgct ctttcgggta gcgcataaaa ccttctgggg gccactgcag
     1621 acagtcatgg aggactacgt gcaacactac gctctagagt gatgatgtgg cagtcaaccg
     1681 agggctgtta agagaagtct gtgcacgaaa atcaatgagt cggtgtaatc ttcaccttcc
     1741 tacaaaccaa ggcttgaggg cctcaatttc gagcaccact tcttgaccaa gccgttcgcc
     1801 ggtataggag ttacgaattc gggtcgtctg ccgtctttga accatagata agtacctccc
     1861 gagtcgggag aacacgaatc gcatcagaac tcatcgaagg gtgttacttt gcccatcggc
     1921 atccccgaag attaattgca ttcgaaattt taacttcgag aatcatccgt gctaaccaac
     1981 gtggacgagc tatgaacgtt tacgtgttca gtatctagga attgacacgt ctccgatagc
     2041 gggagtcacg gcatgagtta cggcttttgt gaccagcgtc tgtgtaattt tattgacagt
     2101 ttcatatata gtgtctcgtc ttgcctttga gattcacagg ctgtaagagg gatgccgtcg
     2161 tgtatatgcc tgcagttaga gcttccttgt caggacgagc aactcactcc ttctggtgat
     2221 gctttgctag catgtagctt cgggtacttg cggaaaagga cagcctacag tgtctgcggc
     2281 ggtacgtcac taatatgaca tacgctgtgg acatcgtgaa gcggcgaatc gctccagatt
     2341 gaggtttaca gagtactcca tagtcacaca actcagaata tgcaggtgct gcagtgaccg
     2401 aatcctctca ggtcgcgctg gcggcccgtt tgagagtgct gaaagcgaag tttatccgca
     2461 ccgtcatcca gtgcacctcc cccagagaat gattacttgt gaatgcggtt aacctgaggg
     2521 ccgacgccgc gattaaattt ttaacgtgat gcacaagctg ggggcactct cccctggggt
     2581 caggttataa ttttaaaccc actacggacc tggactgcaa gatgcttagt ccacgcgtct
     2641 catagaccgt cggtctatag ccatgttcaa cgccgccagg ctgggcgcat aactgtatgt
     2701 gactttgtcc acatccaagg cgacgatggc ggcagttggc cggagggctg gctgtcaagg
     2761 ggccacgcgg atgaagggag ctcatacgca actatggagc gcggtatata gtataggatc
     2821 cgccgggaca cgcgggaatc aagggacaaa ggagcagtca gcgtgaaatc tttacgggtg
     2881 aagcgcgctc gtaccaccta gtcccgccta aggcccttgt gtcgtgcaac atcgctgggg
     2941 cgacggatga caatatgatc aaggcccgag agctttgaat gtaagaggcg tattataaag
     3001 cctaggtccg catgttgata cggggagcag cccatgcgca gtaggaaggt ccgaccaccc
     3061 ttttctgagg tctactccaa cctccttgcc gagtgctctc cgctatctcc acgcatagca
     3121 tactccgctt ctaacatctt cttcaaagca tacacgctac ttagcagagc aaatcgacac
     3181 gcccatggag tggtccgctc aaaccgcgca attaagtata agctaatgtg acggaccgga
     3241 cctacgattt cccccaactc gtggggatca gttgactcgt ctcagcttat gagacataat
     3301 tacggtcata gggtacgcct atgcccctca gagatctcta ggtcatatgc cctatcgggt
     3361 gctcggtgac cattgtctca tacaaagtca tatttggcag ggcttctatg gacgagtatt
     3421 atgtgcccaa aagggaacgt ctgagttggg attcgcttct aaccaatact ctaatggtac
     3481 ctaatttaat gtggtcgctg cgcaccccgc cgggtacaac tgctgagttc ctgagatcgg
     3541 tacggcaaac atcttttatc cagacgaggg ttctggccag aattgatgtt cagatacgag
     3601 aaggcatgtc atgtcacaag gattttggag aatt
''',
        '''LOCUS       part_K                  3634 bp    ds-DNA  circular UNK 28-OCT-2019
DEFINITION  .
ACCESSION   part_K
VERSION     part_K
KEYWORDS    "creator:SynthSys Center" "marker:KanR, BlpR".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     misc_feature    1268..1284
                     /label="feature"
                     /ApEinfo_fwdcolor="#b4abac"
     CDS             complement(join(3199..3634,1..380))
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     misc_feature    2374..2377
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#faac61"
     misc_feature    2395..2419
                     /label="feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     promoter        1382..1560
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     misc_feature    2370..2373
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_feature    2389..2393
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     misc_feature    1239..1261
                     /label="another feature"
                     /ApEinfo_fwdcolor="#f8d3a9"
     rep_origin      complement(2510..3098)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
     terminator      2117..2369
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     CDS             1562..2112
                     /label="another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     rep_origin      671..1106
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#ffef86"
''',
        'Metadata1 for part_K', 'Metadata2 for part_K'
    ),
    (
        'part_L',
        '''ORIGIN
        1 gtccgatgcg agcaccacct atagataccc taatcagcat ggtcccggtg accattacca
       61 ccgacggcct tgttttaaac gaaggcgctc actatacaag aaactaatgt ctcccacaag
      121 gctcagaccc agccggtgct ttcagcgctc cccgttggcg ccccccgaag accatagaga
      181 ggtgctctgt taacaactat aatgaaaggg attattgagg gagtcggagg tgaattctga
      241 agccttaccc ctacggtcgt ggtcatcatc acggatcgcg ataagcggcc ggtctgttac
      301 ggacgccact acgaggaggg tttcttggca atcaaggcta cgcattacaa ccatcggggt
      361 cggtctggct gaacaatctt gatgttcaac ctttgttacc ggctccttgt gttattatcc
      421 ggtttttttg atctatatgt ttataaggaa gcggtgcaga ccgtagagaa ataccggcta
      481 aggcagtcgc tcaatgcaca ccagcatcgc cgagacatat aaagccataa acattaattg
      541 gaatgctagt cacagcgtgt tgactggacg cttaaaagtt tctacaaaag gtatagtttt
      601 aacatttctc tgagttagcg ctatcgggtc gactgacatc tagaattgac acgttactca
      661 cctggatcct agtcacactc cgggaaccgc tcgcttagtg gcagtagcgg ggggcgtccc
      721 taggatagct gcgcaaagcg tggtttatct aatcgtggct tgtgtggaca tgatgctgtt
      781 caggcgcatc tgtggccgcc acagcattac caagctaatc agaaaccgcg gcatgtcccc
      841 gcgatgcgag aggcatctgt tagagccaaa gagtggaggt ctgtaatatc gctagcggag
      901 attttctaag gggcgtggga ctttatcccg atatggctga acgatccaac ggcaaatggt
      961 tctgccattg aactttaacc atgaatagac tccaggcgag ggataatgga aataatagca
     1021 ggaaatgaca atctactagg tgcgcactat caagattgta tgtttgcgta tcgtagttct
     1081 aatggtcttg tgtcccgtat tgggtcggac ggtgctgttc ctgatcgtag caatgcccgg
     1141 ttagaatcac gcatcggagc gatcatttcg ggcgatacgc gttcgttctg tccactcttg
     1201 tgggagtccg acgtatgcct gctcggtatt gttttattca ggcgacagta cctcccattc
     1261 aaatatacag ggttatacta cacaccacga gtaataaaga cttgtgttgt gccaattcct
     1321 cgtgggactg ccacgctcgt gagcttcgcg gttacatacc cctatgtccc caactgttgg
     1381 tggcaaatcc ctaatctaat ctattgtatc acagtcggtc aacgactgaa attggggcaa
     1441 caccaacagt gctcttacga tattcgtctg tgtctaacaa tccaacaaaa agctgttatt
     1501 atgatggagg atcttggatt gcgtcagtac tacggtgttg tacgagaatt gcatgatgat
     1561 tcctcgatgt cgagtccttg caagatcgtg atctatagcc aagctgctgc tgtagatcaa
     1621 cctcaccggg gcgggggtcc cattactatg ggcggcatta gagtaccata gatgaacgcg
     1681 agcagagaca tgtcaatgta tcccagtcct gcagtactac gtcagtctct gctcctgggg
     1741 ccaatgcacg cctgtagaat agggccgtgt tcccatgaga aatttgaggg acttacgatt
     1801 agatggcctt ggcaggggac tccgcgtacg tgtggggatg aaccggcacc ggccttggcc
     1861 ctgagtaaat gtaagcgata tcaacttcgc tattcatgat cagtggtatt cgtgtaggca
     1921 ccgaaagtcc cttggcgaaa gtagagacat acgatacgga ccgcaagcga gcagaatttc
     1981 agcacatgac agttattagc ttttttctaa tggtgctccc attatgattg tcaagtcccg
     2041 actgctagct gtaggatgaa acacatccac acacgtgtta ccatactgat tttgacacta
     2101 cttcctaacc gggccgaatg gtatctgctg cgatgcttag ttacgtaaca gccgaatgtc
     2161 acgccactgc atattgtggt tctggaattc cacctgcata tgtcacgcta gagaccgtat
     2221 gaggtgggcg gatagtggtg acttcaaacc catatgcagg tgctgca
''',
        '''LOCUS       part_L                  2267 bp    ds-DNA  circular UNK 28-OCT-2019
DEFINITION  .
ACCESSION   part_L
VERSION     part_L
KEYWORDS    "creator:SynthSys Center" "marker:KanR".
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     misc_feature    join(2263..2267,1)
                     /label="feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     terminator      2..59
                     /label="another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     CDS             complement(1140..1955)
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#b7e6d7"
     misc_feature    2206..2209
                     /label="feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_feature    2202..2205
                     /label="another feature"
                     /ApEinfo_fwdcolor="#faac61"
     terminator      2139..2182
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
     rep_origin      complement(254..842)
                     /label="feature"
                     /ApEinfo_fwdcolor="#ffef86"
     misc_feature    2248..2251
                     /label="another feature"
                     /ApEinfo_fwdcolor="#ff0000"
     misc_feature    2185..2190
                     /label="yet another feature"
                     /ApEinfo_fwdcolor="#84b0dc"
     terminator      1025..1119
                     /label="feature"
                     /ApEinfo_fwdcolor="#c6c9d1"
''',
        'Metadata1 for part_L', 'Metadata2 for part_L'
    )
])

    conn.commit()
    print("Mock data inserted into 'sample' table.")

    cursor.close()
    conn.close()


if __name__ == "__main__":
    create_db_and_insert_data()
