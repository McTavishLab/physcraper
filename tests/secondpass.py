from physcraper import physcraper_setup, physcraper_scrape
import pickle


runname="asc_test_full"

test3 = physcraper_scrape('{}/{}_scrape.p'.format(runname,runname))
test3.run_blast()
pickle.dump(test3, open('{}/{}_scrape.p'.format(runname,runname), 'wb'))
