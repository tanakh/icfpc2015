import           Control.Applicative
import           Data.List
import           System.Directory
import           System.IO
import           Text.Regex.TDFA

main :: IO ()
main = do
    files <- filter (".json" `isSuffixOf`) <$> getDirectoryContents "out"

    let tt = [ ((ss!!1, ss!!2), (read (ss!!3) :: Int, file))
             | file <- files
             , let (ss:_) = file =~ "^([0-9]+)-([0-9]+)-([0-9]+)\\.json$" :: [[String]]
             ]

        gg =
            sort
            $ map snd
            $ map snd
            $ map maximum
            $ groupBy (\a b -> fst a == fst b)
            $ sort tt

    mapM_ (hPutStrLn stderr) gg

    cons <- mapM (readFile . ("out/" ++)) gg

    putStrLn "["
    putStr $ intercalate "," cons
    putStrLn "]"
